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

# TODO: Do I need to handle QuadraticFunction here?
#function _process_arg(arg::MOI.ScalarQuadraticFunction)::NormalizedExpressionTuple
#    return normalize_function(arg)
#end

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
    #
    # TODO: I think this normalize_function call here is bad for performance.
    # Really, I should preserve small types as much as possible during this
    # recursive function, then normalize once at the very end.
    # I need something like a _recursively_process_nonlinear_function, which
    # is just this function without the call to normalize, then a normalize_function(SNF)
    # that calls the former and normalizes the result.
    #println(processed_args)
    combined_expr = _combine_subexpressions(fcn.head, processed_args)
    #println(combined_expr)
    affine, nonlinear =  normalize_function(combined_expr)
    affine = _combine_like_terms(affine)
    return (affine, nonlinear)
end

const NormalizedNode = Union{
    NormalizedExpressionTuple,
    MOI.ScalarNonlinearFunction, # I think SNF should be allowed here...
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

const VariableAffineNormalizedNode = Union{
    # TODO: Make these types AffineNormalizedNode, and add a
    # ConstantNormalizedNode
    MOI.ScalarAffineFunction,
    MOI.VariableIndex,
}

const NonlinearNormalizedNode = Union{
    NormalizedExpressionTuple,
    MOI.ScalarNonlinearFunction,
    MOI.ScalarQuadraticFunction,
}

#
# Methods to help build up the sub-expression
#
function _get_constant(arg::Real)::Float64
    return arg
end

function _get_constant(arg::MOI.ScalarAffineFunction)::Float64
    return arg.constant
end

function _get_constant(arg::MOI.ScalarQuadraticFunction)::Float64
    return arg.constant
end

function _get_constant(arg::NormalizedExpressionTuple)::Float64
    affine, nonlinear = arg
    return affine.constant
end

function _get_constant(arg::MOI.VariableIndex)::Float64
    return 0.0
end

function _get_linear_terms(arg::Real)::Vector{MOI.ScalarAffineTerm{Float64}}
    return []
end

function _get_linear_terms(arg::MOI.VariableIndex)::Vector{MOI.ScalarAffineTerm{Float64}}
    return [MOI.ScalarAffineTerm(1.0, arg)]
end

function _get_linear_terms(arg::MOI.ScalarAffineFunction)::Vector{MOI.ScalarAffineTerm{Float64}}
    return arg.terms
end

function _get_linear_terms(arg::MOI.ScalarQuadraticFunction)::Vector{MOI.ScalarAffineTerm{Float64}}
    return arg.affine_terms
end

function _get_linear_terms(arg::NormalizedExpressionTuple)::Vector{MOI.ScalarAffineTerm{Float64}}
    affine, nonlinear = arg
    return affine.terms
end

# Note: _get_nonlinear just needs to return a valid ScalarNonlinearFunction node
function _get_nonlinear(arg::Real)
    return nothing
end

function _get_nonlinear(arg::MOI.VariableIndex)
    return nothing
end

function _get_nonlinear(arg::MOI.ScalarAffineFunction)
    return nothing
end

function _get_nonlinear(arg::MOI.ScalarQuadraticFunction)
    new_quadratic = MOI.ScalarQuadraticFunction(
        arg.quadratic_terms,
        Vector{MOI.ScalarAffineTerm{Float64}}(),
        0.0,
    )
    return new_quadratic
end

function _get_nonlinear(arg::NormalizedExpressionTuple)
    affine, nonlinear = arg
    return nonlinear
end

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

function _negate(expr::Real)::Real
    return - expr
end

function _negate(expr::MOI.VariableIndex)::MOI.ScalarAffineFunction
    terms = [MOI.ScalarAffineTerm(- 1.0, expr)]
    affine = MOI.ScalarAffineFunction(terms, 0.0)
    return affine
end

function _negate(expr::MOI.ScalarAffineFunction)::MOI.ScalarAffineFunction
    terms = [
        MOI.ScalarAffineTerm(- term.coefficient, term.variable)
        for term in expr.terms
    ]
    neg_expr = MOI.ScalarAffineFunction(terms, - expr.constant)
    return neg_expr
end

function _negate(expr::MOI.ScalarNonlinearFunction)::MOI.ScalarNonlinearFunction
    return MOI.ScalarNonlinearFunction(:-, [expr])
end

function _negate(expr::NormalizedExpressionTuple)::NormalizedExpressionTuple
    affine, nonlinear = expr
    neg_affine = _negate(affine)
    neg_nonlinear = _negate(nonlinear)
    return (neg_affine, neg_nonlinear)
end

function _combine_normalized_minus_args(args::Vector{Any})::NormalizedNode
    if length(args) == 2
        return _combine_normalized_plus_args(
            Vector{Any}([args[1], _negate(args[2])])
        )
    else
        throw(ArgumentError("- operator must have exactly two arguments"))
    end
end

function _combine_normalized_times_args(args::Vector{Any})::NormalizedNode
    if length(args) < 2
        throw(ArgumentError("* node must have at least two arguments"))
    elseif length(args) > 2
        # n-ary product is tricky to normalize as we have to consider
        # distribution of constant, linear, and nonlinear terms.
        throw(ArgumentError("n-ary product is not supported for now"))
    elseif all(typeof(arg) <: Real for arg in args)
        return prod(args)
    elseif any(typeof(arg) <: Real && arg == 0.0 for arg in args)
        # If any arg is a zero-valued constant, return zero
        return 0.0
    elseif (
        any(typeof(arg) <: NonlinearNormalizedNode for arg in args)
        && any(typeof(arg) <: Real for arg in args)
    )
        # Nonlinear * Real. Just multiply everything by this real
        nl_idx = findfirst(arg -> typeof(arg) <: NonlinearNormalizedNode, args)
        const_idx = findfirst(arg -> typeof(arg) <: Real, args)
        factor = args[const_idx]
        constant = _get_constant(args[nl_idx])
        lin_terms = _get_linear_terms(args[nl_idx])
        nonlinear = _get_nonlinear(args[nl_idx])

        new_constant = factor * constant
        new_terms = MOI.ScalarAffineTerm{Float64}[
            MOI.ScalarAffineTerm(factor * term.coefficient, term.variable)
            for term in lin_terms
        ]
        new_affine = MOI.ScalarAffineFunction(new_terms, constant)
        new_nl = MOI.ScalarNonlinearFunction(:*, [factor, nonlinear])
        return (new_affine, new_nl)
    elseif (
        any(typeof(arg) <: NonlinearNormalizedNode for arg in args)
        && any(typeof(arg) <: VariableAffineNormalizedNode for arg in args)
    )
        # Nonlinear * Affine.
        nl_idx = findfirst(arg -> typeof(arg) <: NonlinearNormalizedNode, args)
        aff_idx = findfirst(arg -> typeof(arg) <: VariableAffineNormalizedNode, args)

        # Note that either constant can be zero
        nl_const = _get_constant(args[nl_idx])
        aff_const = _get_constant(args[aff_idx])
        new_const = nl_const * aff_const

        nl_lin_terms = _get_linear_terms(args[nl_idx])
        aff_lin_terms = _get_linear_terms(args[aff_idx])
        # NOTE: This would be a lot easier with AffineExpression and operator overloading.
        new_terms = Vector{MOI.ScalarAffineTerm{Float64}}()
        for term in nl_lin_terms
            if term.coefficient != 0.0 && aff_const != 0.0
                coef = term.coefficient * aff_const
                new_term = MOI.ScalarAffineTerm(coef, term.variable)
                push!(new_terms, new_term)
            end
        end
        for term in aff_lin_terms
            if term.coefficient != 0.0 && nl_const != 0.0
                coef = term.coefficient * nl_const
                new_term = MOI.ScalarAffineTerm(coef, term.variable)
                push!(new_terms, new_term)
            end
        end
        new_affine = MOI.ScalarAffineFunction(new_terms, new_const)

        nl_nonlinear = _get_nonlinear(args[nl_idx])
        nonlinear_terms = Vector{Any}()
        if aff_const != 0.0
            nl = MOI.ScalarNonlinearFunction(:*, [aff_const, nl_nonlinear])
            push!(nonlinear_terms, nl)
        end
        # Assume we have at least one linear term in the affine function
        affine_lin_function = MOI.ScalarAffineFunction(aff_lin_terms, 0.0)
        if length(nl_lin_terms) != 0
            nl_lin_function = MOI.ScalarAffineFunction(nl_lin_terms, 0.0)
            nl = MOI.ScalarNonlinearFunction(:*, [nl_lin_function, affine_lin_function])
            push!(nonlinear_terms, nl)
        end
        nl = MOI.ScalarNonlinearFunction(:*, [nl_nonlinear, affine_lin_function])
        push!(nonlinear_terms, nl)
        new_nonlinear = MOI.ScalarNonlinearFunction(:+, nonlinear_terms)
        return (new_affine, new_nonlinear)
    elseif all(typeof(arg) <: NonlinearNormalizedNode for arg in args)
        # Nonlinear * nonlinear
        constants = [_get_constant(arg) for arg in args]
        lin_terms = [_get_linear_terms(arg) for arg in args]
        nl_fcns = [_get_nonlinear(arg) for arg in args]

        new_constant = constants[1] * constants[2]

        new_lin_terms = Vector{MOI.ScalarAffineTerm{Float64}}()
        if length(lin_terms[1]) != 0 && constants[2] != 0.0
            for term in lin_terms[1]
                coef = constants[2] * term.coefficient
                new_term = MOI.ScalarAffineTerm(coef, term.variable)
                push!(new_lin_terms, new_term)
            end
        end
        if length(lin_terms[2]) != 0 && constants[1] != 0.0
            for term in lin_terms[2]
                coef = constants[1] * term.coefficient
                new_term = MOI.ScalarAffineTerm(coef, term.variable)
                push!(new_lin_terms, new_term)
            end
        end
        new_affine = MOI.ScalarAffineFunction(new_lin_terms, new_constant)

        new_nl_terms = Vector{Any}()
        if constants[1] != 0.0
            nl = MOI.ScalarNonlinearFunction(:*, [constants[1], nl_fcns[2]])
            push!(new_nl_terms, nl)
        end
        if constants[2] != 0.0
            nl = MOI.ScalarNonlinearFunction(:*, [constants[2], nl_fcns[1]])
            push!(new_nl_terms, nl)
        end
        if length(lin_terms[1]) != 0
            lin = MOI.ScalarAffineFunction(lin_terms[1], 0.0)
            nl = MOI.ScalarNonlinearFunction(:*, [lin, nl_fcns[2]])
            push!(new_nl_terms, nl)
        end
        if length(lin_terms[2]) != 0
            lin = MOI.ScalarAffineFunction(lin_terms[2], 0.0)
            nl = MOI.ScalarNonlinearFunction(:*, [lin, nl_fcns[1]])
            push!(new_nl_terms, nl)
        end
        if length(lin_terms[1]) != 0 && length(lin_terms[2]) != 0
            lin1 = MOI.ScalarAffineFunction(lin_terms[1], 0.0)
            lin2 = MOI.ScalarAffineFunction(lin_terms[2], 0.0)
            nl = MOI.ScalarNonlinearFunction(:*, [lin1, lin2])
            push!(new_nl_terms, nl)
        end
        nl = MOI.ScalarNonlinearFunction(:*, [nl_fcns[1], nl_fcns[2]])
        push!(new_nl_terms, nl)
        new_nonlinear = MOI.ScalarNonlinearFunction(:+, new_nl_terms)
        return (new_affine, new_nonlinear)
    # FIXME: This does not reduce properly.
    #elseif any(typeof(arg) <: NonlinearNormalizedNode for arg in args)
    #    # Arg could be NormalizedExpressionTuple or ScalarQuadraticFunction
    #    #
    #    # If we have any nonlinear expression (that is not cancelled out by
    #    # a zero), we return a nonlinear expression.
    #    args = filter(arg -> !(typeof(arg) <: Real && arg == 1.0), args)
    #    constants = [_get_constant(arg) for arg in args]
    #    lin_terms = [_get_linear_terms(arg) for arg in args]
    #    nl_terms = [_get_nonlinear(arg) for arg in args]
    #    # NOTE: Here is where we assume that we only have two arguments.
    #    # Otherwise we would need more complicated distribution.
    #    constant = prod(constants; init = 1.0)
    #    terms = Vector{MOI.ScalarAffineTerm{Float64}}()
    #    for term in lin_terms[1]
    #        new_term = MOI.ScalarAffineTerm{Float64}(constants[2] * term.coefficient)
    #        push!(terms, new_term)
    #    end
    #    for term in lin_terms[2]
    #        new_term = MOI.ScalarAffineTerm{Float64}(constants[1] * term.coefficient)
    #        push!(terms, new_term)
    #    end
    #    # TODO: Combine nonlinear terms
    #    #
    #    #affine, _ = normalize_function(0.0)
    #    #nonlinear = MOI.ScalarNonlinearFunction(:*, args)
    #    #return (affine, nonlinear)
    elseif count(arg -> (typeof(arg) <: VariableAffineNormalizedNode), args, init=0) == 1
        # We only have one linear/affine subexpression
        aff_idx = findfirst(arg -> (typeof(arg) <: VariableAffineNormalizedNode), args)
        affine = args[aff_idx]
        factor = prod(
            arg for (idx, arg) in enumerate(args) if (idx != aff_idx && arg != 1.0);
            init = 1.0,
        )
        constant = 0.0
        terms = MOI.ScalarAffineTerm{Float64}[]
        constant += _collect_affine!(terms, affine)
        # Note: on rare occasions this may simplify further into a simple variable,
        # in the case `affine` is just a variable and `factor == 1`
        new_terms = [
            MOI.ScalarAffineTerm(factor * term.coefficient, term.variable)
            for term in terms
        ]
        new_const = factor * constant
        return MOI.ScalarAffineFunction(new_terms, new_const)
        # Under what conditions do we return affine?
        # Note that input args have been normalized.
        # Have affine expressions necessarily been maximally reduced? With my current
        # code, I don't think so. I would like to rely on this, however, so that affine
        # times affine is always nonlinear.
    # TODO: Still need to implement affine * affine
    else
        throw(ArgumentError("Unhandled * argument"))
    end
end

_OPERATOR_DISPATCHER = Dict{Symbol, Function}(
    :+ => _combine_normalized_plus_args,
    :- => _combine_normalized_minus_args,
    :* => _combine_normalized_times_args,
)

function _combine_subexpressions(
    op::Symbol, args::Vector{Any}
)::NormalizedNode
    return _OPERATOR_DISPATCHER[op](args)
end
