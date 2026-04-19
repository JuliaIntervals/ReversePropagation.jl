# Interface for using ChainRules.jl

# Base.conj(x::Num) = x   # assuming reals
# Base.complex(x::Num) = x
# Base.float(x::Num) = x

@scalar_rule(^(x::Num, n::Integer), (n==1 ? 1 : n==2 ? 2x : n*x^(n-1), ZeroTangent()))

# The default ChainRules rule for abs uses sign(x), which isn't reversible
# in this package's op set. Replace with x/abs(x), whose ops (/, abs) both
# have reverse contractors.
@scalar_rule(abs(x::Num), x / abs(x))

# Default ChainRules rules for exp2/exp10/log2/log10 reference
# `IrrationalConstants.Logtwo`/`Logten`, which IntervalArithmetic only
# supports for `Base.MathConstants` irrationals. Replace with plain
# Float64 log values so the generated SSA is interval-friendly.
@scalar_rule(exp2(x::Num),  exp2(x)  * log(2.0))
@scalar_rule(exp10(x::Num), exp10(x) * log(10.0))
@scalar_rule(log2(x::Num),  inv(x * log(2.0)))
@scalar_rule(log10(x::Num), inv(x * log(10.0)))

adj(f, z̄::Num, x::Num) = (rrule(f, x)[2](z̄))[2]
adj(f, z̄, x) = adj(f, Num(z̄), Num(x))

adj(f, z̄::Num, x, y) = simplify( (rrule(f, x, y)[2])(z̄)[2:end] )
adj(f, z̄, x, y) = adj(f, Num(z̄), Num(x), Num(y))

struct Tangent{T <: Real} <: Real
end

# tangent(var) = Variable(Symbol(value(var), "̇"))
# tangent(var::Num) = Sym{Tangent{Real}}(Symbol(var, "̇"))
tangent(var::Num) = tangent(value(var))
tangent(var::Sym) = Sym{Tangent{Real}}(Symbol(var, "̇"))
tangent(x::Real) = 0  # derivative of a constant

function tangent(eq::Assignment)

    vars = args(eq)
    rhs_tangents = Num.(tangent.(args(eq)))
    lhs_tangent = Num(tangent(lhs(eq)))

    return Assignment(lhs_tangent,
                    frule( (ZeroTangent(), rhs_tangents...),
                    op(eq), vars...)[2]
    )
end



# julia> frule((ZeroTangent(), ẋ, ẏ), *, x, y)
# (x * y, (ẋ * y) + (x * ẏ))

# julia> dargs = Num.(dotify.(args(eq)))

# julia> a = Num.(args(eq))
# 2-element Vector{Num}:
#  x
#  y

# julia> frule( (ZeroTangent(), dargs...), op(eq), a...)
# (x * y, (ẋ * y) + (x * ẏ))

