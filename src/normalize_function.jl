#  ___________________________________________________________________________
#
#  MathProgIncidence.jl: Math Programming Incidence Graph Analysis
#  Copyright (c) 2023. Triad National Security, LLC. All rights reserved.
#
#  This program was produced under U.S. Government contract 89233218CNA000001
#  for Los Alamos National Laboratory (LANL), which is operated by Triad
#  National Security, LLC for the U.S. Department of Energy/National Nuclear
#  Security Administration. All rights in the program are reserved by Triad
#  National Security, LLC, and the U.S. Department of Energy/National Nuclear
#  Security Administration. The Government is granted for itself and others
#  acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license
#  in this material to reproduce, prepare derivative works, distribute copies
#  to the public, perform publicly and display publicly, and to permit others
#  to do so.
#
#  This software is distributed under the 3-clause BSD license.
#  ___________________________________________________________________________

"""
Methods to normalize an MOI scalar function

"""

import JuMP
import MathOptInterface as MOI

const NormalizedExpressionTuple = Tuple{
    MOI.ScalarAffineFunction,
    MOI.ScalarNonlinearFunction,
}

struct NormalizedExpression
    affine::Union{MOI.ScalarAffineFunction, Nothing}
    nonlinear::Union{MOI.ScalarNonlinearFunction, Nothing}
end

function normalize_function(fcn::MOI.ScalarAffineFunction)::NormalizedExpressionTuple
    return (fcn, MOI.ScalarNonlinearFunction(:+, [0.0]))
end

function normalize_function(fcn::MOI.ScalarQuadraticFunction)::NormalizedExpressionTuple
    affine = MOI.ScalarAffineFunction(fcn.affine_terms, fcn.constant)
    # TODO: should I check the quadratic_term type and make sure to use it for the
    # empty affine subexpression and the constant term? How do I do this?
    quadratic = MOI.ScalarQuadraticFunction(
        fcn.quadratic_terms,
        Vector{MOI.ScalarAffineTerm{Float64}}(),
        0.0,
    )
    nonlinear = MOI.ScalarNonlinearFunction(:+, [quadratic])
    return (affine, nonlinear)
end

function normalize_function(fcn::MOI.ScalarNonlinearFunction)::NormalizedExpressionTuple
    processed_args = Vector{Any}()
    for arg in fcn.args
        if _is_leaf(arg)
            # Note that we still may want to process ScalarQuadraticFunction and
            # ScalarAffineFunction to detect whether they resolve to a simpler
            # form.
            push!(processed_args, arg)
        else
            # arg <: ScalarNonlinearFunction
            # We need to construct a new normalized expression here. Without
            # considering all arguments simultaneously, we can't determine
            # how to build up the current affine and nonlinear subexpressions.
            # And a normalized expression is the data structure that lets us
            # process the arguments.
            # TODO: Would it be easier if we convert all arguments into
            # normalized expressions? Rather than just ScalarNonlinearFunction?
            #
            # Imagine *(0, SNF). This can be handled more efficiently with zero
            # a constant. We don't have to check whether its normalized expr
            # resolves to a trivial value.
            #
            # should normalize_function be able to return a simpler type (e.g.
            # constant) to help with the case where the SNF reduces?
            # The user-facing function should be type-stable, implying that this
            # recursive procesing should use its own function?
            norm_subexpr = normalize_function(arg)
            push!(processed_args, norm_subexpr)
        end
    end
    # Now I just have to implement _combine_subexpressions! for every operator
    return _combine_subexpressions!(fcn.head, processed_args)
end
