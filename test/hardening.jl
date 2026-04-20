using IntervalArithmetic
using ReversePropagation: binary_functions, unary_functions,
                          cse_equations, gradient_code, binarize_ssa,
                          _rev_binary_lookup, _rev_unary_lookup,
                          DIFF_HELPERS

# Ops that may appear in a gradient SSA: functions with reverse contractors
# (registry) plus the hand-written derivative helpers (`_ln2`, `_abs_subgrad`
# etc.) that appear on rule rhs's.
const _registered_ops = union(
    keys(_rev_binary_lookup),
    keys(_rev_unary_lookup),
    DIFF_HELPERS,
)

function _is_reversible(eq)
    rhs = Symbolics.value(eq.rhs)
    SymbolicUtils.iscall(rhs) || return true
    op = SymbolicUtils.operation(rhs)
    return op in _registered_ops
end

@testset "derivative hardening" begin

    @testset "registry shape" begin
        # `max` / `min` are binary; keep them out of the unary list.
        @test !(:max in unary_functions)
        @test !(:min in unary_functions)
        @test :max in keys(binary_functions)
        @test :min in keys(binary_functions)
        # `sign` is fully excluded: distributional derivative at 0, so it
        # is not an allowed op in user expressions.
        @test !(:sign in unary_functions)
        @test_throws ErrorException ReversePropagation.unary_rule(sign, 0.0)
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
        # For every function in the current registries *that has a scalar
        # rule*, verify that the binarized gradient SSA contains only ops
        # that are either in the reversibility registry or are declared
        # `DIFF_HELPERS`.
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

    @testset "abs subgradient over intervals" begin
        @variables x
        g = ReversePropagation.gradient(abs(x), [x])
        # Real eval: sign(x).
        @test g((3.0,))[2] == (1.0,)
        @test g((-2.5,))[2] == (-1.0,)
        @test g((0.0,))[2] == (0.0,)
    end
end
