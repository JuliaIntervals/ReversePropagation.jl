"""
    SSAFunction

A symbolic function in SSA (single static assignment) form.

Fields:
- `code`: vector of `Assignment`s, each of the form `_a := f(x, y)`
- `output`: the symbolic variable holding the final result

An `SSAFunction` is produced by `cse_equations` (common subexpression elimination)
and can be passed to `forward_backward_code` or `gradient_code` as input.
"""
struct SSAFunction
    code::Vector{Assignment}
    output
end

function Base.show(io::IO, ssa::SSAFunction)
    n = length(ssa.code)
    print(io, "SSAFunction($n assignments, output = $(ssa.output))")
end

function Base.show(io::IO, ::MIME"text/plain", ssa::SSAFunction)
    n = length(ssa.code)
    println(io, "SSAFunction with $n assignments:")
    for eq in ssa.code
        println(io, "  $(eq.lhs) := $(eq.rhs)")
    end
    print(io, "Output: $(ssa.output)")
end
