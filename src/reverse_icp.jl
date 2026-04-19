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
                    :max   => :max_rev,
                    :min   => :min_rev,
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


# `sign` is kept in the reversibility registry (IntervalContractors has
# `sign_rev`), so downstream interval-reverse consumers can propagate
# through it, but it is deliberately *not* given a scalar rule in
# `scalar_rules.jl` — its derivative is distributional and should not be
# taken. It reaches gradient SSA only via `_abs_subgrad`'s runtime
# dispatch for Real inputs.
const unary_functions = [:sqrt, :abs,
            :exp, :exp2, :exp10, :expm1,
            :log, :log2, :log10, :log1p,
            :sin, :cos, :tan,
            :asin, :acos, :atan,
            :sinh, :cosh, :tanh,
            :asinh, :acosh, :atanh,
            :inv, :sign];

const _rev_unary_lookup = Dict{Function,Function}()

for f in unary_functions
    f_rev = Symbol(f, :_rev)
    rev_sym = Symbol(:_rev_, f)
    @eval $rev_sym(z, x) = $f_rev(z, x)
    @eval @register_symbolic $rev_sym(z, x) false
    _rev_unary_lookup[eval(f)] = eval(rev_sym)
end

rev(f_val, z, x) = _rev_unary_lookup[f_val](z, x)


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

    reverse_code = rev.(reverse(forward_code), Ref(params))

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

# ex = x^2 + a * y^2
# C = forward_backward_contractor(ex, vars, [a])

# const CC2 = forward_backward_contractor(ex, vars, [a])

# CC2((-10..10, -10..10), 0..1, 5)
# @btime CC2((-10..10, -10..10), 0..1, 6)
# @btime CC2((-10..10, -10..10), 0..1, 7)


# using BenchmarkTools

# @btime CC((-10..10, -10..10), 0..1)

# @code_native C((-10..10, -10..10), 0..1)

# CC(IntervalBox(-10..10, 2), 0..1)



# vars = @variables x, y
# ex = 3x^2 + 4x^2 * y


# code, final, return_tuple = forward_backward_code(vars, 3x^2 + 4x^2 * y)

# expr = forward_backward_expr(vars, ex)

# C = forward_backward_contractor(vars, ex)


# using IntervalArithmetic

# C( (1..2, 3..4) )

# ex = x^2 + y^2

# C( (-10..10, -10..10) )
