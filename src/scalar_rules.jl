"""
Hand-written scalar derivative rules for the functions in
`binary_functions` / `unary_functions`.

Why not ChainRules or DiffRules:

- Our registered op set is small (~22 unary, 5 binary). A local table is
  easier to audit than filtering upstream rules.
- Every rhs in this table is constrained to use only ops that themselves
  have interval reverse contractors in the registry, so the binarized
  output of `gradient_code` stays reversible end-to-end.
- `abs`, `max`, `min` need interval-subgradient reverses eventually;
  extending this table is a localised change, whereas overriding
  ChainRules/DiffRules entries is monkey-patching.
"""

# Partials ∂f/∂x for each unary function in the registry.
const UNARY_RULES = Dict{Function, Function}(
    sqrt  => x -> inv(2 * sqrt(x)),
    abs   => x -> x / abs(x),              # reversible surrogate for sign(x)
    exp   => x -> exp(x),
    exp2  => x -> exp2(x)  * log(2.0),     # avoid IrrationalConstants.Logtwo
    exp10 => x -> exp10(x) * log(10.0),
    expm1 => x -> exp(x),
    log   => x -> inv(x),
    log2  => x -> inv(x * log(2.0)),
    log10 => x -> inv(x * log(10.0)),
    log1p => x -> inv(1 + x),
    sin   => x -> cos(x),
    cos   => x -> -sin(x),
    tan   => x -> inv(cos(x)^2),           # sec²(x)
    asin  => x -> inv(sqrt(1 - x^2)),
    acos  => x -> -inv(sqrt(1 - x^2)),
    atan  => x -> inv(1 + x^2),
    sinh  => x -> cosh(x),
    cosh  => x -> sinh(x),
    tanh  => x -> 1 - tanh(x)^2,
    asinh => x -> inv(sqrt(x^2 + 1)),
    acosh => x -> inv(sqrt(x^2 - 1)),
    atanh => x -> inv(1 - x^2),
    inv   => x -> -inv(x)^2,
)

# Partials (∂f/∂x, ∂f/∂y) for each binary function in the registry.
const BINARY_RULES = Dict{Function, Function}(
    (+) => (x, y) -> (1, 1),
    (-) => (x, y) -> (1, -1),
    (*) => (x, y) -> (y, x),
    (/) => (x, y) -> (inv(y), -x * inv(y)^2),
    (^) => (x, y) -> (y * x^(y - 1), x^y * log(x)),
)

function unary_rule(f, x)
    rule = get(UNARY_RULES, f, nothing)
    isnothing(rule) &&
        error("no scalar derivative rule registered for unary $(f)")
    return rule(x)
end

function binary_rule(f, x, y)
    rule = get(BINARY_RULES, f, nothing)
    isnothing(rule) &&
        error("no scalar derivative rule registered for binary $(f)")
    return rule(x, y)
end

# Reverse-mode adjoint: given the cotangent `z̄` of the output of `f`, return
# the cotangent contribution(s) to the input(s).
adj(f, z̄::Num, x::Num) = unary_rule(f, x) * z̄
adj(f, z̄, x)           = adj(f, Num(z̄), Num(x))

function adj(f, z̄::Num, x::Num, y::Num)
    dx, dy = binary_rule(f, x, y)
    return (dx * z̄, dy * z̄)
end
adj(f, z̄, x, y) = adj(f, Num(z̄), Num(x), Num(y))
