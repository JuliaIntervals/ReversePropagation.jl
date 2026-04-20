@testset "forward-mode tangent" begin

    vars = @variables x, y

    @testset "directional derivative matches reverse-mode gradient" begin
        # f = x² + y² at (3, 4): ∇f = (2x, 2y) = (6, 8)
        f = x^2 + y^2
        t = tangent(f, [x, y])
        g = ReversePropagation.gradient(f, [x, y])

        pt = (3.0, 4.0)
        (v_fwd_x, dx) = t(pt, (1.0, 0.0))
        (v_fwd_y, dy) = t(pt, (0.0, 1.0))
        (v_rev, grad) = g(pt)

        @test v_fwd_x == v_rev == 25.0
        @test v_fwd_y == v_rev
        @test (dx, dy) == grad
    end

    @testset "transcendental" begin
        # d/dx sin(x*y) = y cos(x*y) ; d/dy sin(x*y) = x cos(x*y)
        t = tangent(sin(x*y), [x, y])
        pt = (1.0, 2.0)
        (_, dx) = t(pt, (1.0, 0.0))
        (_, dy) = t(pt, (0.0, 1.0))
        @test dx ≈ 2.0 * cos(2.0)
        @test dy ≈ 1.0 * cos(2.0)
    end

    @testset "arbitrary direction" begin
        # Full Gateaux derivative: directional derivative along v equals ⟨∇f, v⟩
        f = x^2 * y + y^3
        t = tangent(f, [x, y])
        g = ReversePropagation.gradient(f, [x, y])

        pt = (2.0, 3.0)
        dir = (0.7, -0.4)
        (_, dfdv) = t(pt, dir)
        (_, grad) = g(pt)
        @test dfdv ≈ grad[1] * dir[1] + grad[2] * dir[2]
    end

    @testset "tangent_code surfaces input/output tangents" begin
        ssa = ReversePropagation.tangent_code(x^2 + y^2, [x, y])
        @test length(ssa.variables.input_tangents) == 2
        @test ssa.variables.output_tangent isa Num
        # Every original forward assignment produces exactly one tangent assignment.
        fwd_ssa = ReversePropagation.cse_equations(x^2 + y^2)
        @test length(ssa.code) == 2 * length(fwd_ssa.code)
    end
end
