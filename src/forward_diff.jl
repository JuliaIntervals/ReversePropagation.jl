"""
Forward-mode automatic differentiation on top of `SSAFunction`.

Given an SSA `_a := f(x₁, …, xₙ)`, the tangent (directional derivative)
is ``_ȧ = Σᵢ (∂f/∂xᵢ) · ẋᵢ``. We reuse the scalar rule table from
`scalar_rules.jl` — the same rules power both reverse- and forward-mode,
so there is no second source of truth to keep in sync.
"""


"Symbolic tangent variable for `v` (overdot notation)."
dotted(v::Num) = dotted(Symbolics.value(v))
dotted(v::Real) = Num(0)
function dotted(v::SymbolicUtils.BasicSymbolic)
    return Num(Symbolics.variable(Symbol(string(v), '̇')))
end


"Tangent assignment for one forward assignment `eq`."
function tangent_eq(eq::Assignment)
    vs = args(eq)
    lhs_dot = dotted(lhs(eq))
    if length(vs) == 1
        dfdx = unary_rule(op(eq), vs[1])
        rhs = dfdx * dotted(vs[1])
    elseif length(vs) == 2
        x, y = vs
        dfdx, dfdy = binary_rule(op(eq), x, y)
        rhs = dfdx * dotted(x) + dfdy * dotted(y)
    else
        error("tangent_eq: unsupported arity $(length(vs))")
    end
    return Assignment(lhs_dot, rhs)
end


"""
    tangent_code(ssa::SSAFunction, vars) -> SSAFunction
    tangent_code(ex::Num, vars) -> SSAFunction

Extend `ssa` (or `cse_equations(ex)`) with forward-mode tangent assignments.

The returned `SSAFunction.code` contains the original forward assignments
followed by one tangent assignment per forward assignment. The extra
outputs are surfaced in `variables`:

- `input_tangents`: the tangent inputs the caller must supply (one per `vars`).
- `output_tangent`: the tangent of `ssa.output`.
"""
function tangent_code(ssa::SSAFunction, vars)
    fwd = ssa.code
    tang = [tangent_eq(eq) for eq in fwd]
    combined = [fwd; tang]
    extras = (
        input_tangents = collect(dotted.(vars)),
        output_tangent = dotted(ssa.output),
    )
    return SSAFunction(combined, ssa.output, merge(ssa.variables, extras))
end

tangent_code(ex::Num, vars) = tangent_code(cse_equations(ex), vars)


"""
    tangent(ex, vars) -> callable

Build a compiled forward-mode function of the form
`(__inputs, __tangents) -> (value, directional_derivative)`.

The directional derivative is ``Σᵢ (∂f/∂xᵢ)(__inputs) · __tangents[i]`` —
set `__tangents = (1, 0, …, 0)` for ``∂f/∂x₁``, and so on.
"""
function tangent(ex::Num, vars)
    ssa = tangent_code(ex, vars)

    code = Expr(:block, toexpr.(ssa.code)...)

    input_vars         = toexpr(Symbolics.MakeTuple(vars))
    input_tangent_vars = toexpr(Symbolics.MakeTuple(ssa.variables.input_tangents))
    final              = toexpr(ssa.output)
    final_tangent      = toexpr(ssa.variables.output_tangent)

    # NB: the first arg is named `__inputs` (not `__args`); `__args` collides
    # with RuntimeGeneratedFunctions' own internal binding and triggers a
    # spurious `BoundsError` at call time.
    full_code = :(
        (__inputs, __tangents) -> begin
            $input_vars         = __inputs
            $input_tangent_vars = __tangents
            $code
            return ($final, $final_tangent)
        end
    )

    return @RuntimeGeneratedFunction(full_code)
end

tangent(f, vars) = tangent(f(vars), vars)
