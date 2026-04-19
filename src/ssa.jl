"""
    SSAFunction

A symbolic function in SSA (single static assignment) form.

Fields:
- `code`: vector of `Assignment`s, each of the form `_a := f(x, y)`
- `output`: the symbolic variable holding the final result
- `variables`: a `NamedTuple` of extra named variables (e.g. `constraint`, `gradient`)

An `SSAFunction` is produced by `cse_equations` (common subexpression elimination)
and can be passed to `forward_backward_code` or `gradient_code` as input.
"""
struct SSAFunction
    code::Vector{Assignment}
    output::Num
    variables::NamedTuple
end

SSAFunction(code, output) = SSAFunction(code, output, (;))

function Base.show(io::IO, ssa::SSAFunction)
    n = length(ssa.code)
    print(io, "SSAFunction($n assignments, output = $(ssa.output))")
end

_show_lhs(x) = string(x)
_show_lhs(t::Symbolics.MakeTuple) = "(" * join(t.elems, ", ") * ")"

function _show_rhs(rhs)
    val = Symbolics.value(rhs)
    if SymbolicUtils.iscall(val)
        f = SymbolicUtils.operation(val)
        name = string(nameof(f))
        if startswith(name, "_rev_")
            op = name[6:end]
            args = join(SymbolicUtils.arguments(val), ", ")
            return "rev($op)($args)"
        end
    end
    return string(rhs)
end

function Base.show(io::IO, ::MIME"text/plain", ssa::SSAFunction)
    n = length(ssa.code)
    println(io, "SSAFunction with $n assignments:")
    for eq in ssa.code
        println(io, "  $(_show_lhs(eq.lhs)) := $(_show_rhs(eq.rhs))")
    end
    print(io, "Output: $(ssa.output)")
    if !isempty(keys(ssa.variables))
        for (k, v) in pairs(ssa.variables)
            print(io, "\n", k, ": ", v)
        end
    end
end


"""
    binarize_ssa(ssa::SSAFunction) -> SSAFunction

Decompose any n-ary or compound-argument assignments in `ssa.code` into
elementary unary/binary operations on atomic arguments.

Useful as a post-processing step on the output of `gradient_code`: the
reverse-mode adjoint pass can produce assignments like
`_ā := _b̄ * cos(_a)` or `_t := *(2, x, _ā)`, which downstream consumers
(e.g. an interval reverse-propagation pipeline) need split into
elementary ops — `_tmp := cos(_a); _ā := _b̄ * _tmp` etc.

Preserves `ssa.output` and `ssa.variables`; only rewrites `ssa.code`.
"""
function binarize_ssa(ssa::SSAFunction)
    binarized = _binarize_ssa(ssa.code)
    return SSAFunction(binarized, ssa.output, ssa.variables)
end

function _binarize_ssa(code)
    result = Assignment[]
    for eq in code
        rhs_val = Symbolics.value(eq.rhs)
        if SymbolicUtils.iscall(rhs_val) && _needs_decomposition(rhs_val)
            sub_ssa = cse_equations(eq.rhs)
            append!(result, sub_ssa.code[1:end-1])
            push!(result, Assignment(eq.lhs, sub_ssa.code[end].rhs))
        else
            push!(result, eq)
        end
    end
    return result
end

function _needs_decomposition(ex)
    args = SymbolicUtils.arguments(ex)
    length(args) > 2 && return true
    return any(SymbolicUtils.iscall, args)
end
