
lhs(eq::Assignment) = eq.lhs
rhs(eq::Assignment) = value(eq.rhs)

op(eq::Assignment) = operation(rhs(eq))


_map_to_num(x) = Num(x)
_map_to_num(x::ExactReal) = x

# args(eq::Assignment) = Num.(arguments(rhs(eq)))
args(eq::Assignment) = _map_to_num.(arguments(rhs(eq)))

name(var) = value(var).name
bar(var) = variable(Symbol(var, '̄'))  # the character is the overbar symbol (on top of first `'`)




# "Construct adjoint of a given Assignment with accumulation"
# function adj(eq::Assignment)
#     vars = args(eq)
#     adjoints = Adjoint(op(eq), vars)

#     bar_lhs = bar(lhs(eq))

#     eqns = map(zip(vars, adjoints)) do (var, adj)
#         barred = bar(var)

#         Assignment(barred, barred + bar_lhs * adj)
#     end

#     return eqns
# end


"Make a new variable by adding a number to an old variable, e.g. `x` -> `x3`"
function add_number(var, num)
    return variable(var, num)   # Variable(Symbol(var, "_", num))
end


function numbered_variable!(num_times_used, var, increment=false)
    is_new = false

    if !haskey(num_times_used, var)
        is_new = true

        num = 0
        num_times_used[var] = 0
    else

        if increment
            num_times_used[var] += 1
        end

        num = num_times_used[var]
    end

    return is_new, add_number(var, num)
end


function adj(eq::Assignment)
    vars = value.(args(eq))
    # vars = value.(_map_to_num.(arguments(eq)))

    bar_lhs = bar(lhs(eq))

    adjoints = adj(op(eq), bar_lhs, vars...)

    eqns = map(zip(vars, adjoints)) do (var, adj)
        #if var isa Sym{Tangent{Real}}  # for linearization pass
            barred = bar(var)

            Assignment(barred, barred + adj)
        #end
    end

    return eqns
end

# modifies the Set assigned of variables which have already been assigned
function simple_adj(eq::Assignment, assigned)

    vars = ReversePropagation.args(eq)
    bar_lhs = bar(lhs(eq))

    adjoints = adj(op(eq), bar_lhs, vars...)

    eqns = Assignment[]

    for (var, adjoint) in zip(vars, adjoints)

        value = Symbolics.value(var)
        ((value isa Real) || (value isa Sym)) && continue

        barred = bar(var)

        if any(x -> isequal(x, barred), assigned)
            push!(eqns, Assignment(barred, barred + adjoint))

        else
            push!(eqns, Assignment(barred, adjoint))
            push!(assigned, barred)
        end

    end

    return eqns
end


function adj(num_times_used, eq::Assignment)

    # @show num_times_used
    # @show eq

    vars = args(eq)
    # adjoints = Adjoint(op(eq), vars)

    # is_new, bar_lhs = numbered_variable!(num_times_used, bar(lhs(eq)))

    bar_lhs = bar(lhs(eq))
    is_new, numbered_bar_lhs = numbered_variable!(num_times_used, bar_lhs)

    adjoints = adj(op(eq), numbered_bar_lhs, vars...)

    eqns = Assignment[]

    # @show bar_lhs
    for (var, adj) in zip(vars, adjoints)

        # @show typeof(var)
        !isa(var, Num) && continue

        barred = bar(var)

        is_new, numbered = numbered_variable!(num_times_used, barred)

        if is_new # first use
            # Assignment(numbered, numbered_bar_lhs * adj)
            push!(eqns, Assignment(numbered, adj))

        else  # re-use so generate new variable
            new_num, barred_new = numbered_variable!(num_times_used, barred, true)
            # Assignment(barred_new, numbered + numbered_bar_lhs * adj)
            push!(eqns, Assignment(barred_new, numbered + adj))

        end
    end

    return eqns
end


function reverse_pass(vars, code, final)
    num_times_used = Dict{Num, Int}()

    num, final_bar = numbered_variable!(num_times_used, bar(final))
    reverse_code = [Assignment(final_bar, 1)]

    for eq in reverse(code)
        # @show eq
        append!(reverse_code, adj(num_times_used, eq))
    end

    final_vars = [numbered_variable!(num_times_used, bar(var))[2] for var in vars]
    return reverse_code, final_vars
end

function simple_reverse_pass(vars, forward_code)

    assigned = Set()  # which variables have already been assigned

    reverse_code =  reduce(vcat, simple_adj.(reverse(forward_code), Ref(assigned)))

    final_vars = bar.(vars)

    return reverse_code, final_vars, assigned

end



"""
Return code for forward and reverse pass, as a vector of Assignments.
`final` is the output variable from the forward pass.
`gradient_vars` are the output variables from the reverse pass.
"""
function gradient_code(ex, vars)
    ssa = cse_equations(ex)
    return gradient_code(ssa, vars)
end

function gradient_code(ssa::SSAFunction, vars)
    forward_code = ssa.code
    final = ssa.output

    reverse_code, gradient_vars, assigned = simple_reverse_pass(vars, forward_code)

    initialization_code = [Assignment(bar(final), 1)]

    unassigned = setdiff(gradient_vars, assigned)
    append!(initialization_code, [Assignment(var, 0) for var in unassigned])

    code = [forward_code; initialization_code; reverse_code]

    return SSAFunction(code, final, (; gradient=gradient_vars))
end




make_tuple(args) = Expr(:tuple, args...)
make_tuple(s::Symbol) = make_tuple([s])

# toexpr(ex::Assignment) = toexpr(Equation(ex.lhs, ex.rhs))

function gradient_expr(ex, vars)
    result_ssa = gradient_code(ex, vars)

    code = Expr(:block, toexpr.(result_ssa.code)...)

    return code, result_ssa.output, result_ssa.variables.gradient
end


function gradient(ex::Num, vars)

    code, final_var, gradient_vars = gradient_expr(ex, vars)

    input_vars = toexpr(Symbolics.MakeTuple(vars))
    final = toexpr(final_var)
    gradient = toexpr(Symbolics.MakeTuple(gradient_vars))

    full_code = :(
        (__args,) -> begin
            $input_vars = __args
            $code
            return $(final), $(gradient)
        end
    )

    return @RuntimeGeneratedFunction(full_code)

end



gradient(f, vars) = gradient(f(vars), vars)

