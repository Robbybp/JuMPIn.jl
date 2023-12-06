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

# TODO: Should normalize_function be in charge of removing duplicates, i.e. linear
# cancellations in affine expressions? Probably, for now, there should be one
# function for normalizing and one function for consolidating terms.
#
# I think we actually need to greedily consolidate terms. Consider
# (x - x) * fcn(y)
# This doesn't simplify to 0 if we don't consolidate (x - x) greedily.
# Does this matter? This will potentially cause incorrect incidence in the
# nonlinear subexpression, but we don't guarantee that nonlinear incidence
# will correctly deal with cancellations.
# Is there a case where lack of cancellation leads to incorrect linear incidence?

function _combine_like_terms(
    affine::MOI.ScalarAffineFunction
)::MOI.ScalarAffineFunction
    # TODO: Use "eltype" of ScalarAffineFunction here
    var_to_coef = Dict{MOI.VariableIndex, Float64}()
    for term in affine.terms
        var = term.variable
        if var in keys(var_to_coef)
            var_to_coef[var] += term.coefficient
        else
            var_to_coef[var] = term.coefficient
        end
    end
    terms = [MOI.ScalarAffineTerm(coef, var) for (var, coef) in var_to_coef]
    new_affine = MOI.ScalarAffineFunction(terms, affine.constant)
    return new_affine
end

function normalize_function(fcn::MOI.VariableIndex)::NormalizedExpressionTuple
    terms = [MOI.ScalarAffineTerm(1.0, fcn)]
    affine = MOI.ScalarAffineFunction(terms, 0.0)
    nonlinear = MOI.ScalarNonlinearFunction(:+, [0.0])
    return (affine, nonlinear)
end

function normalize_function(fcn::Real)::NormalizedExpressionTuple
    terms = Vector{MOI.ScalarAffineTerm{Float64}}()
    affine = MOI.ScalarAffineFunction(terms, fcn)
    nonlinear = MOI.ScalarNonlinearFunction(:+, [0.0])
    return (affine, nonlinear)
end

function normalize_function(fcn::MOI.ScalarAffineFunction)::NormalizedExpressionTuple
    return (
        _combine_like_terms(fcn),
        MOI.ScalarNonlinearFunction(:+, [0.0]),
    )
end

function normalize_function(fcn::MOI.ScalarQuadraticFunction)::NormalizedExpressionTuple
    affine = MOI.ScalarAffineFunction(fcn.affine_terms, fcn.constant)
    affine = _combine_like_terms(affine)
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

function normalize_function(fcn::NormalizedExpressionTuple)::NormalizedExpressionTuple
    return fcn
end

function _process_arg(arg)
    return arg
end

function _process_arg(arg::MOI.ScalarNonlinearFunction)::NormalizedExpressionTuple
    return normalize_function(arg)
end

function normalize_function(fcn::MOI.ScalarNonlinearFunction)::NormalizedExpressionTuple
    processed_args = Vector{Any}()
    for arg in fcn.args
        push!(processed_args, _process_arg(arg))
        #if _is_leaf(arg)
        #    # Note that we still may want to process ScalarQuadraticFunction and
        #    # ScalarAffineFunction to detect whether they resolve to a simpler
        #    # form.
        #    push!(processed_args, arg)
        #else
        #    # arg <: ScalarNonlinearFunction
        #    # We need to construct a new normalized expression here. Without
        #    # considering all arguments simultaneously, we can't determine
        #    # how to build up the current affine and nonlinear subexpressions.
        #    # And a normalized expression is the data structure that lets us
        #    # process the arguments.
        #    # TODO: Would it be easier if we convert all arguments into
        #    # normalized expressions? Rather than just ScalarNonlinearFunction?
        #    #
        #    # Imagine *(0, SNF). This can be handled more efficiently with zero
        #    # a constant. We don't have to check whether its normalized expr
        #    # resolves to a trivial value.
        #    #
        #    # should normalize_function be able to return a simpler type (e.g.
        #    # constant) to help with the case where the SNF reduces?
        #    # The user-facing function should be type-stable, implying that this
        #    # recursive procesing should use its own function?
        #    norm_subexpr = normalize_function(arg)
        #    push!(processed_args, norm_subexpr)
        #end
    end
    # Now I just have to implement _combine_subexpressions! for every operator
    affine, nonlinear =  normalize_function(_combine_subexpressions(fcn.head, processed_args))
    affine = _combine_like_terms(affine)
    return (affine, nonlinear)
end

const NormalizedNode = Union{
    NormalizedExpressionTuple,
    MOI.ScalarAffineFunction,
    MOI.ScalarQuadraticFunction,
    MOI.VariableIndex,
    Real,
}

const AffineNormalizedNode = Union{
    MOI.ScalarAffineFunction,
    MOI.VariableIndex,
    Real,
}

const NonlinearNormalizedNode = Union{
    NormalizedExpressionTuple,
    MOI.ScalarQuadraticFunction,
}

#
# Methods to help build up the affine sub-expression
#
function _collect_affine!(
    terms::Vector{MOI.ScalarAffineTerm{Float64}},
    arg::Real,
)::Float64
    return arg
end

function _collect_affine!(
    terms::Vector{MOI.ScalarAffineTerm{Float64}},
    arg::MOI.VariableIndex,
)::Float64
    constant = 0.0
    push!(terms, MOI.ScalarAffineTerm(1.0, arg))
    return constant
end

function _collect_affine!(
    terms::Vector{MOI.ScalarAffineTerm{Float64}},
    arg::MOI.ScalarAffineFunction,
)::Float64
    for term in arg.terms
        push!(terms, term)
    end
    return arg.constant
end

function _collect_affine!(
    terms::Vector{MOI.ScalarAffineTerm{Float64}},
    arg::MOI.ScalarQuadraticFunction,
)::Float64
    for term in arg.affine_terms
        push!(terms, term)
    end
    return arg.constant
end

function _collect_affine!(
    terms::Vector{MOI.ScalarAffineTerm{Float64}},
    arg::NormalizedExpressionTuple,
)::Float64
    affine, nonlinear = arg
    return _collect_affine!(terms, affine)
end

# TODO: _collect_nonlinear! methods

function _collect_nonlinear!(subexpr::Vector{Any}, arg::Real)
end

function _collect_nonlinear!(subexpr::Vector{Any}, arg::MOI.VariableIndex)
end

function _collect_nonlinear!(subexpr::Vector{Any}, arg::MOI.ScalarAffineFunction)
end

function _collect_nonlinear!(subexpr::Vector{Any}, arg::MOI.ScalarQuadraticFunction)
    quadratic = MOI.ScalarQuadraticFunction(
        arg.quadratic_terms,
        Vector{MOI.ScalarAffineTerm{Float64}}(),
        0.0,
    )
    push!(subexpr, quadratic)
end

function _collect_nonlinear!(subexpr::Vector{Any}, arg::NormalizedExpressionTuple)
    affine, nonlinear = arg
    push!(subexpr, nonlinear)
end

#function _combine_normalized_plus_args(args::Vector{NormalizedNode})::NormalizedNode
function _combine_normalized_plus_args(args::Vector{Any})::NormalizedNode
    # This is actualy more complicated if we expect all args to be normalized expressions
    # What is the simplest version of this function I could write? One that operates
    # on some NormalizedNode, which is a union of SNF, SAF, SQF, VarIndex, and Real.
    # Then it also returns a NormalizedNode?
    if length(args) == 0
        return 0.0
    elseif length(args) == 1
        return args[1]
    elseif all(typeof(arg) <: Real for arg in args)
        return sum(args)
    elseif all(typeof(arg) <: AffineNormalizedNode for arg in args)
        terms = Vector{MOI.ScalarAffineTerm{Float64}}()
        constant = 0.0
        for arg in args
            # Note that we are building up terms in-place
            constant += _collect_affine!(terms, arg)
        end
        return MOI.ScalarAffineFunction(terms, constant)
    elseif any(typeof(arg) <: NonlinearNormalizedNode for arg in args)
        terms = Vector{MOI.ScalarAffineTerm{Float64}}()
        constant = 0.0
        # Args in nonlinear_subexpr can be ScalarNonlinearFunction or
        # ScalarQuadraticFunction
        nonlinear_subexpr = Vector{Any}()
        for arg in args
            constant += _collect_affine!(terms, arg)
            _collect_nonlinear!(nonlinear_subexpr, arg)
        end
        affine = MOI.ScalarAffineFunction(terms, constant)
        nonlinear = MOI.ScalarNonlinearFunction(:+, nonlinear_subexpr)
        return (affine, nonlinear)
    end
end

function _combine_normalized_minus_args(args::Vector{Any})::NormalizedNode
    # Note that this function is inefficient as it creates an entire
    # NormalizedExpressionTuple even if the args just contains a single variable.
    affine, nonlinear = normalize_function(_combine_normalized_plus_args(args))
    affine_terms = [
        MOI.ScalarAffineTerm(- term.coefficient, term.variable)
        for term in affine.terms
    ]
    new_affine = MOI.ScalarAffineFunction(affine_terms, - affine.constant)
    new_nonlinear = MOI.ScalarNonlinearFunction(:-, nonlinear.args)
    return (new_affine, new_nonlinear)
end

_OPERATOR_DISPATCHER = Dict{Symbol, Function}(
    :+ => _combine_normalized_plus_args,
    :- => _combine_normalized_minus_args,
)

function _combine_subexpressions(
    op::Symbol, args::Vector{Any}
)::NormalizedNode
    return _OPERATOR_DISPATCHER[op](args)
end
