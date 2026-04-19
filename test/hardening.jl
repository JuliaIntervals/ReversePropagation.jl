using IntervalArithmetic
using ReversePropagation: binary_functions, unary_functions,
                          cse_equations, gradient_code, binarize_ssa,
                          _rev_binary_lookup, _rev_unary_lookup

# The set of operations that can appear in an SSA and have a matching
# reverse contractor (so downstream interval-reverse-propagation consumers
# can rev them without surprises).
const _registered_ops = Set{Function}(
    reduce(vcat, [
        collect(keys(_rev_binary_lookup)),
        collect(keys(_rev_unary_lookup)),
    ])
)

# An assignment is "reversible" if its rhs is a constant, a bare variable,
# or a call whose operation is in `_registered_ops`.
function _is_reversible(eq)
    rhs = Symbolics.value(eq.rhs)
    SymbolicUtils.iscall(rhs) || return true
    op = SymbolicUtils.operation(rhs)
    return op in _registered_ops
end

@testset "derivative hardening" begin

    @testset "registry shape" begin
        # `sign` is distributional at 0 and has no `sign_rev`; omitting it
        # is load-bearing for the smooth-ops guarantee below.
        @test !(:sign in unary_functions)
        # `max` / `min` are binary and require subgradient reverses.
        @test !(:max in unary_functions)
        @test !(:min in unary_functions)
    end

    @testset "binarize_ssa" begin
        @variables x y
        # `3x + x*y` has an n-ary `*` (or `+`) after CSE normalization
        # that should survive binarization as binary ops only.
        ssa = cse_equations(3x + x*y)
        bin = binarize_ssa(ssa)
        @test bin isa SSAFunction
        @test bin.output === ssa.output
        # Every call-rhs in the binarized SSA has at most two atomic args.
        for eq in bin.code
            rhs = Symbolics.value(eq.rhs)
            SymbolicUtils.iscall(rhs) || continue
            args = SymbolicUtils.arguments(rhs)
            @test length(args) <= 2
            @test all(a -> !SymbolicUtils.iscall(Symbolics.value(a)), args)
        end
    end

    @testset "gradient SSA uses only registered ops" begin
        # For every function in the current registries, verify that the
        # binarized gradient SSA contains only ops that have a reverse
        # contractor. This is the contract that lets downstream
        # derivative-based contractors consume `gradient_code` output
        # without hitting `sign_rev not defined` etc.
        @variables x y

        for f_sym in unary_functions
            f = eval(f_sym)
            ssa = cse_equations(f(x))
            grad = binarize_ssa(gradient_code(ssa, [x]))
            @test all(_is_reversible, grad.code)
        end

        for (f_sym, _) in binary_functions
            f = eval(f_sym)
            ssa = cse_equations(f(x, y))
            grad = binarize_ssa(gradient_code(ssa, [x, y]))
            @test all(_is_reversible, grad.code)
        end
    end
end
