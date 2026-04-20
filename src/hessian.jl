"""
Hessian computation via forward-over-reverse: build the gradient SSA
(one reverse-mode pass), then attach a tangent (forward-mode) pass that
propagates derivatives of every SSA variable — including the gradient
components `variables.gradient`. Evaluating the tangent pass with the
`i`-th unit vector as the input direction yields the `i`-th row of the
Hessian.
"""


"""
    hessian(ex, vars) -> callable

Build a compiled function `__inputs -> (value, gradient, Hessian)` where
`gradient` is an `n`-tuple and `Hessian` is an `n × n` matrix.

Implementation: one compiled `(inputs, direction) -> (value, gradient,
gradient_tangent)` built by forward-over-reverse, invoked `n` times at
runtime with unit-vector directions. Cost is `n` forward/reverse passes
for the Hessian build, which is acceptable for the small-`n` symbolic
regime this package targets; denser cases can compile `n` tangent passes
into a single SSA later.
"""
function hessian(ex::Num, vars)
    n = length(vars)

    grad_ssa = gradient_code(cse_equations(ex), vars)
    g_syms   = grad_ssa.variables.gradient

    # Binarize before tangenting: the adjoint pass can emit n-ary muls
    # (e.g. `2 * _b̄ * x`) that `tangent_eq` doesn't handle. `binarize_ssa`
    # preserves `.variables` (so `g_syms` still point into the code), it
    # just rewrites compound rhs's into elementary binary/unary form.
    grad_ssa = binarize_ssa(grad_ssa)

    # Attach tangent pass. This produces tangent assignments for every SSA
    # variable, including the `g_syms`. The tangent of `g_syms[j]` under
    # input direction `eᵢ` is `H[i, j]`.
    tang_ssa = tangent_code(grad_ssa, vars)
    grad_tangent_syms = [dotted(g) for g in g_syms]

    code                = Expr(:block, toexpr.(tang_ssa.code)...)
    input_vars          = toexpr(Symbolics.MakeTuple(vars))
    input_tangent_vars  = toexpr(Symbolics.MakeTuple(tang_ssa.variables.input_tangents))
    value_sym           = toexpr(tang_ssa.output)
    grad_tuple          = toexpr(Symbolics.MakeTuple(g_syms))
    grad_tangent_tuple  = toexpr(Symbolics.MakeTuple(grad_tangent_syms))

    inner = @RuntimeGeneratedFunction(:(
        (__inputs, __tangents) -> begin
            $input_vars         = __inputs
            $input_tangent_vars = __tangents
            $code
            return ($value_sym, $grad_tuple, $grad_tangent_tuple)
        end
    ))

    return function hessian_closure(__inputs)
        T = typeof(first(__inputs))
        one_T  = one(T)
        zero_T = zero(T)

        H  = Matrix{T}(undef, n, n)
        local value
        local grad

        for i in 1:n
            dir = ntuple(k -> k == i ? one_T : zero_T, n)
            v, g, gt = inner(__inputs, dir)
            if i == 1
                value = v
                grad  = g
            end
            for j in 1:n
                H[i, j] = gt[j]
            end
        end

        return (value, grad, H)
    end
end

hessian(f, vars) = hessian(f(vars), vars)
