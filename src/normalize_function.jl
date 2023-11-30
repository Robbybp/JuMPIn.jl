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
    # Here I need a function that recursively builds affine and nonlinear subexpressions
    #
    # Again, it seems wrong to enforce the type of these terms' coefficients
    affine_terms = Vector{MOI.ScalarAffineTerm{Float64}}()

    # Initialize with trivial affine and nonlinear expressions
    affine = MOI.ScalarAffineFunction(affine_terms, 0.0)
    nonlinear = MOI.ScalarNonlinearFunction(:+, [0.0])
    _collect_subexpressions!(fcn, affine, nonlinear)
    return (affine, nonlinear)
end

# _collect_subexpressions! needs to be implemented for every node that can appear
# in ScalarNonlinearFunction.
#
# The following won't work, as it doesn't know whether this affine function is
# being being multiplied by some variable expression.
# I need to handle every combination of operator and argument type.
# I will likely need to build up expression from leaf-to-root, although I may
# be able to do some top-down processing for efficiency.
function _collect_subexpressions!(
    fcn::MOI.ScalarAffineFunction,
    affine::MOI.ScalarAffineFunction,
    nonlinear::MOI.ScalarQuadraticFunction,
)
    for term in fcn.terms
        add_to_expression!(affine, term)
    end
    add_to_expression!(affine, fcn.constant)
    return (affine, nonlinear)
end

function _collect_subexpressions!(
    fcn::MOI.ScalarQuadraticFunction,
    affine::MOI.ScalarAffineFunction,
    nonlinear::MOI.ScalarQuadraticFunction,
)
end
