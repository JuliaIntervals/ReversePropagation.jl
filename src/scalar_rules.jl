"""
Hand-written scalar derivative rules for the functions in
`binary_functions` / `unary_functions`.

Each rule is a method of `unary_rule` / `binary_rule`, dispatched on the
function value (`::typeof(sin)` etc.) rather than looked up in a `Dict`.
This is both cheaper (no function-identity hashing at SSA generation) and
more idiomatic — users can extend the rule table by adding their own
methods rather than mutating a module-level dict.

The rules for non-smooth functions (`abs`, `max`, `min`) evaluate to
interval subgradients when the runtime input is an interval, and to
ordinary derivatives for plain `Real`. The rules for
`exp2` / `exp10` / `log2` / `log10` use a type-preserving `_ln2` / `_ln10`
helper so `log(2)` and `log(10)` are computed rigorously in the interval
flavour of the input rather than snapped to a `Float64` literal.
"""


# ---------------------------------------------------------------------------
# Runtime helpers used inside the scalar rules
# ---------------------------------------------------------------------------

# Type-preserving `log(n)` so derivative evaluations stay rigorous over
# intervals. Plain `log(2.0)` would snap to `Float64`, which is fine for
# scalar gradients but lossy for `Interval{T}` / `BareInterval{T}` inputs.
_ln_const(::Type{T}, n::Integer) where {T <: Real} = log(T(n))
_ln_const(::Type{Interval{T}},     n::Integer) where {T} = log(interval(T, n))
_ln_const(::Type{BareInterval{T}}, n::Integer) where {T} = log(bareinterval(T, n))

_ln2(x)  = _ln_const(typeof(x), 2)
_ln10(x) = _ln_const(typeof(x), 10)

@register_symbolic _ln2(x)  false
@register_symbolic _ln10(x) false


# Subgradient of `|x|`:
#   Real:      sign(x)     (0 at the origin)
#   Interval:  {+1}/{-1} if strictly signed, else [-1, 1]
_abs_subgrad(x::Real) = sign(x)

function _abs_subgrad(x::Interval{T}) where {T}
    inf(x) > 0 && return interval(T, 1)
    sup(x) < 0 && return interval(T, -1)
    return interval(T, -1, 1)
end

function _abs_subgrad(x::BareInterval{T}) where {T}
    inf(x) > 0 && return bareinterval(T, 1)
    sup(x) < 0 && return bareinterval(T, -1)
    return bareinterval(T, -1, 1)
end

@register_symbolic _abs_subgrad(x) false


# Subgradients of `max(x, y)` / `min(x, y)` w.r.t. each argument. For Real
# we pick 1 on the ≥/≤ branch (so `max(x, x)` credits x). For Interval, if
# the two ranges are separated we return a sharp 0/1; if they overlap we
# return the full subgradient interval `[0, 1]`.
_max_subgrad_x(x::Real, y::Real) = x >= y ? one(x) : zero(x)
_max_subgrad_y(x::Real, y::Real) = y >  x ? one(y) : zero(y)
_min_subgrad_x(x::Real, y::Real) = x <= y ? one(x) : zero(x)
_min_subgrad_y(x::Real, y::Real) = y <  x ? one(y) : zero(y)

function _max_subgrad_x(x::Interval{T}, y::Interval{T}) where {T}
    inf(x) > sup(y) && return interval(T, 1)
    sup(x) < inf(y) && return interval(T, 0)
    return interval(T, 0, 1)
end
function _max_subgrad_y(x::Interval{T}, y::Interval{T}) where {T}
    inf(y) > sup(x) && return interval(T, 1)
    sup(y) < inf(x) && return interval(T, 0)
    return interval(T, 0, 1)
end
function _min_subgrad_x(x::Interval{T}, y::Interval{T}) where {T}
    sup(x) < inf(y) && return interval(T, 1)
    inf(x) > sup(y) && return interval(T, 0)
    return interval(T, 0, 1)
end
function _min_subgrad_y(x::Interval{T}, y::Interval{T}) where {T}
    sup(y) < inf(x) && return interval(T, 1)
    inf(y) > sup(x) && return interval(T, 0)
    return interval(T, 0, 1)
end

@register_symbolic _max_subgrad_x(x, y) false
@register_symbolic _max_subgrad_y(x, y) false
@register_symbolic _min_subgrad_x(x, y) false
@register_symbolic _min_subgrad_y(x, y) false


# The set of helpers above behave like ops in SSA (they appear on the rhs
# of a gradient assignment) but aren't user-facing math operations. Record
# them so the hardening test can whitelist them when checking that every
# rhs is "reversible".
const DIFF_HELPERS = Set{Function}([
    _ln2, _ln10,
    _abs_subgrad,
    _max_subgrad_x, _max_subgrad_y,
    _min_subgrad_x, _min_subgrad_y,
])


# ---------------------------------------------------------------------------
# Scalar derivative rules, dispatched on function value
# ---------------------------------------------------------------------------

unary_rule(::typeof(sqrt),  x) = inv(2 * sqrt(x))
unary_rule(::typeof(abs),   x) = _abs_subgrad(x)
unary_rule(::typeof(exp),   x) = exp(x)
unary_rule(::typeof(exp2),  x) = exp2(x)  * _ln2(x)
unary_rule(::typeof(exp10), x) = exp10(x) * _ln10(x)
unary_rule(::typeof(expm1), x) = exp(x)
unary_rule(::typeof(log),   x) = inv(x)
unary_rule(::typeof(log2),  x) = inv(x * _ln2(x))
unary_rule(::typeof(log10), x) = inv(x * _ln10(x))
unary_rule(::typeof(log1p), x) = inv(1 + x)
unary_rule(::typeof(sin),   x) = cos(x)
unary_rule(::typeof(cos),   x) = -sin(x)
unary_rule(::typeof(tan),   x) = inv(cos(x)^2)
unary_rule(::typeof(asin),  x) = inv(sqrt(1 - x^2))
unary_rule(::typeof(acos),  x) = -inv(sqrt(1 - x^2))
unary_rule(::typeof(atan),  x) = inv(1 + x^2)
unary_rule(::typeof(sinh),  x) = cosh(x)
unary_rule(::typeof(cosh),  x) = sinh(x)
unary_rule(::typeof(tanh),  x) = 1 - tanh(x)^2
unary_rule(::typeof(asinh), x) = inv(sqrt(x^2 + 1))
unary_rule(::typeof(acosh), x) = inv(sqrt(x^2 - 1))
unary_rule(::typeof(atanh), x) = inv(1 - x^2)
unary_rule(::typeof(inv),   x) = -inv(x)^2

unary_rule(f, x) =
    error("no scalar derivative rule registered for unary $(f)")


binary_rule(::typeof(+),   x, y) = (1, 1)
binary_rule(::typeof(-),   x, y) = (1, -1)
binary_rule(::typeof(*),   x, y) = (y, x)
binary_rule(::typeof(/),   x, y) = (inv(y), -x * inv(y)^2)
binary_rule(::typeof(^),   x, y) = (y * x^(y - 1), x^y * log(x))
binary_rule(::typeof(max), x, y) = (_max_subgrad_x(x, y), _max_subgrad_y(x, y))
binary_rule(::typeof(min), x, y) = (_min_subgrad_x(x, y), _min_subgrad_y(x, y))

binary_rule(f, x, y) =
    error("no scalar derivative rule registered for binary $(f)")


# ---------------------------------------------------------------------------
# Reverse-mode adjoint via the rules
# ---------------------------------------------------------------------------

adj(f, z̄::Num, x::Num) = unary_rule(f, x) * z̄
adj(f, z̄, x)           = adj(f, Num(z̄), Num(x))

function adj(f, z̄::Num, x::Num, y::Num)
    dx, dy = binary_rule(f, x, y)
    return (dx * z̄, dy * z̄)
end
adj(f, z̄, x, y) = adj(f, Num(z̄), Num(x), Num(y))
