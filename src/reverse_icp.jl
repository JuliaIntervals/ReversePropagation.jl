# Reverse interval constraint propagation using Symbolics


# import Base: ⊓
# import Base: ⊔


#

⊔
⊔

# @register a ⊓ b
# @register a ⊔ b

# Possibly should replace all calls to `rev` with calls to the actual
# reverse functions instead for speed


function remove_constant(s)
    value = Symbolics.value(s)
    return (value isa Real || value isa Sym) ? variable(:_) : value
end

remove_parameters(s, params) = (any(x -> isequal(x, s), params)) ? variable(:_) : s


function rev(eq::Assignment, params)

    vars = tuple(args(eq)...)

    # vars = tuple(arguments(eq)...)
    # vars = tuple(_map_to_num.(arguments(eq))...)
    return_vars = remove_constant.(tuple(lhs(eq), vars...))
    # return_vars = remove_constant.(tuple(lhs(eq), vars...))

    return_vars = remove_parameters.(return_vars, Ref(params))

    reverse = rev(op(eq), lhs(eq), vars...)

    return Assignment(Symbolics.MakeTuple(return_vars), reverse)

end


# difference between reverse mode AD and reverse propagation:
# reverse mode AD introduces *new* variables
# reverse propagation can use the *same* variables

# x ~ x ⊓ (z - y)


# reverse ops from IntervalContractors:


const binary_functions = Dict(
                    :+     => :plus_rev,
                    :-     => :minus_rev,
                    :*     => :mul_rev,
                    :/     => :div_rev,
                    :^     => :power_rev,
                    );

# Create individually-named symbolic functions to avoid method overwriting
# (@register_symbolic widens typeof(f) to a single symbolic type, so all
# binary/unary rev methods would collide if registered under the same name)
const _rev_binary_lookup = Dict{Function,Function}()

for (f, f_rev) in binary_functions
    rev_sym = Symbol(:_rev_, f)
    @eval $rev_sym(z, x, y) = $f_rev(z, x, y)
    @eval @register_symbolic $rev_sym(z, x, y) false
    _rev_binary_lookup[eval(f)] = eval(rev_sym)
end

rev(f_val, z, x, y) = _rev_binary_lookup[f_val](z, x, y)


const unary_functions = [:sqrt, :abs,
            :exp, :exp2, :exp10, :expm1,
            :log, :log2, :log10, :log1p,
            :sin, :cos, :tan,
            :asin, :acos, :atan,
            :sinh, :cosh, :tanh,
            :asinh, :acosh, :atanh,
            :inv];

const _rev_unary_lookup = Dict{Function,Function}()

for f in unary_functions
    f_rev = Symbol(f, :_rev)
    rev_sym = Symbol(:_rev_, f)
    @eval $rev_sym(z, x) = $f_rev(z, x)
    @eval @register_symbolic $rev_sym(z, x) false
    _rev_unary_lookup[eval(f)] = eval(rev_sym)
end

rev(f_val, z, x) = _rev_unary_lookup[f_val](z, x)


# When the gradient SSA is routed through `forward_backward_contractor`
# (e.g. `d/dx(x+y) = 1`), the output variable can be a plain `Real` seed
# left by `gradient_code` that then flows into `last ⊓ constraint_var`.
# Fill in the missing methods so the intersect promotes the scalar to the
# same interval flavour as the constraint.
IntervalArithmetic.intersect_interval(x::Real, y::Interval{T}) where {T} =
    intersect_interval(interval(T, x), y)
IntervalArithmetic.intersect_interval(x::Interval{T}, y::Real) where {T} =
    intersect_interval(x, interval(T, y))
IntervalArithmetic.intersect_interval(x::Real, y::BareInterval{T}) where {T} =
    intersect_interval(bareinterval(T, x), y)
IntervalArithmetic.intersect_interval(x::BareInterval{T}, y::Real) where {T} =
    intersect_interval(x, bareinterval(T, y))


"Generate code (as Symbolics.Assignment) for forward--backward (HC4Revise) contractor from a symbolic expression"
function forward_backward_code(ex, vars, params=[])
    ssa = cse_equations(ex)
    return forward_backward_code(ssa, vars, params)
end

"Generate forward--backward contractor code from an SSAFunction"
function forward_backward_code(ssa::SSAFunction, vars, params=[])

    final_var = make_variable(:value)  # to record the output of running the forward interval function
    constraint_var = make_variable(:constraint)   # symbolic constraint variable

    forward_code = ssa.code
    last = ssa.output

    constraint_code = [Assignment(final_var, last),
                        Assignment(last, last ⊓ constraint_var)]

    # Only reverse assignments that are function calls;
    # constant/identity assignments (e.g. from gradient_code initialization)
    # are needed in the forward pass but cannot be reversed.
    reversible_code = filter(eq -> iscall(value(eq.rhs)), forward_code)
    reverse_code = rev.(reverse(reversible_code), Ref(params))

    code = [forward_code; constraint_code; reverse_code]

    return SSAFunction(code, final_var, (; constraint=constraint_var))
end

# code, final, constraint = forward_backward_code(x^2 + a*y^2, [x, y], [a])

# code

Symbolics.toexpr(t::Tuple) = Symbolics.toexpr(Symbolics.MakeTuple(t))

# vars = @variables x, y

# ex = x^2 + a*y^2
# code, final, constraint = forward_backward_code(ex, vars, [a])

# code
# final
# constraint

# dump(toexpr.(code))

# function forward_backward_expr(vars, ex)
#     symbolic_code, final, constraint = forward_backward_code(vars, ex)

#     code = toexpr.(symbolic_code)

#     return_tuple = toexpr(vars)
#     # push!(code.args, :(return $return_tuple))

#     return code, final, return_tuple
# end

"Build Julia code for forward_backward contractor from a symbolic expression"
function forward_backward_expr(ex, vars, params=[])
    ssa = cse_equations(ex)
    return forward_backward_expr(ssa, vars, params)
end

"Build Julia code for forward_backward contractor from an SSAFunction"
function forward_backward_expr(ssa::SSAFunction, vars, params=[])
    result_ssa = forward_backward_code(ssa, vars, params)

    code = toexpr.(result_ssa.code)
    all_code = Expr(:block, code...)

    return all_code, result_ssa.output, result_ssa.variables.constraint
end


function forward_backward_contractor(ssa::SSAFunction, vars, params=[])
    code, final_var, constraint_var = forward_backward_expr(ssa, vars, params)

    input_vars = toexpr(Symbolics.MakeTuple(vars))
    final = toexpr(final_var)
    constraint = toexpr(constraint_var)

    if !isempty(params)
        params_tuple = toexpr(Symbolics.MakeTuple(params))

        function_code = :(
            (__inputs, __constraint, __params) -> begin
                $input_vars = __inputs
                $constraint = __constraint
                $params_tuple = __params
                $code
                return $input_vars, $(final)
            end
        )
    else
        function_code = :(
            (__inputs, __constraint) -> begin
                $input_vars = __inputs
                $constraint = __constraint
                $code
                return $input_vars, $(final)
            end
        )
    end

    return @RuntimeGeneratedFunction(function_code)
end

function forward_backward_contractor(ex, vars, params=[])
    ssa = cse_equations(ex)
    return forward_backward_contractor(ssa, vars, params)
end


"""
    derivative_ssa(ex, diff_var, vars)

Build an `SSAFunction` that computes the partial derivative ∂ex/∂diff_var
using reverse-mode AD (gradient_code).

The returned SSA can be passed to `forward_backward_contractor` to create
a contractor for a derivative constraint.
"""
function derivative_ssa(ex, diff_var, vars)
    ssa = cse_equations(ex)
    grad_ssa = gradient_code(ssa, vars)

    idx = findfirst(v -> isequal(v, diff_var), vars)
    isnothing(idx) && error("Variable $diff_var not found in vars $vars")

    # The adjoint pass may produce n-ary operations (e.g. *(2, x, _b̄)).
    # Decompose them into binary ops so rev() can handle them.
    binarized = _binarize_ssa(grad_ssa.code)

    # Promote the integer seed from `gradient_code` (`_ā := 1`, any
    # `unassigned := 0`) to intervals. Without this, identity-propagating
    # gradients such as `d/dx(x-y) = 1` leave behind `Int`-valued entries
    # that cause the reverse pass (e.g. `mul_rev`) to recurse forever via
    # its `mul_rev(a,b,c) = mul_rev(promote(a,b,c)...)` fallback.
    promoted = _promote_integer_seeds(binarized)

    return SSAFunction(promoted, grad_ssa.variables.gradient[idx])
end

function _promote_integer_seeds(code)
    result = Assignment[]
    for eq in code
        rhs_val = Symbolics.value(eq.rhs)
        if rhs_val isa Integer
            push!(result, Assignment(eq.lhs, interval(rhs_val)))
        else
            push!(result, eq)
        end
    end
    return result
end

"""
Ensure all SSA assignments are fully decomposed into elementary (unary/binary)
operations on atomic arguments (variables and constants, no compound subexpressions).
This is needed because gradient_code may produce assignments like `_ā := _b̄*cos(_a)`
which must be split into `_tmp := cos(_a); _ā := _b̄*_tmp` for `rev()` to work.
"""
function _binarize_ssa(code)
    result = Assignment[]
    for eq in code
        rhs_val = value(eq.rhs)
        if iscall(rhs_val) && _needs_decomposition(rhs_val)
            # Decompose via CSE into elementary assignments
            sub_ssa = cse_equations(eq.rhs)
            # Add intermediate assignments, and redirect the last one
            # to the original LHS (avoids creating an identity copy)
            append!(result, sub_ssa.code[1:end-1])
            push!(result, Assignment(eq.lhs, sub_ssa.code[end].rhs))
        else
            push!(result, eq)
        end
    end
    return result
end

"Check if an expression has >2 arguments or any compound (non-atomic) arguments"
function _needs_decomposition(ex)
    args = arguments(ex)
    length(args) > 2 && return true
    return any(iscall, args)
end

"""
    derivative_contractor(ex, diff_var, vars, params=[])

Build a forward--backward contractor for the partial derivative ∂ex/∂diff_var.
"""
function derivative_contractor(ex, diff_var, vars, params=[])
    ssa = derivative_ssa(ex, diff_var, vars)
    return forward_backward_contractor(ssa, vars, params)
end
