@testset "hessian" begin

    vars = @variables x y z

    @testset "quadratic: x² + y²" begin
        # ∇f = (2x, 2y), Hessian = [[2, 0], [0, 2]]
        H = hessian(x^2 + y^2, [x, y])
        val, grad, M = H((3.0, 4.0))
        @test val == 25.0
        @test grad == (6.0, 8.0)
        @test M == [2.0 0.0; 0.0 2.0]
    end

    @testset "cross term: x * y" begin
        # ∇f = (y, x), Hessian = [[0, 1], [1, 0]]
        H = hessian(x * y, [x, y])
        _, _, M = H((2.0, 3.0))
        @test M == [0.0 1.0; 1.0 0.0]
    end

    @testset "polynomial: x²y + y³" begin
        # ∇f = (2xy, x² + 3y²)
        # Hessian:
        #   ∂²f/∂x² = 2y
        #   ∂²f/∂x∂y = 2x
        #   ∂²f/∂y² = 6y
        H = hessian(x^2 * y + y^3, [x, y])
        val, grad, M = H((2.0, 3.0))
        @test val ≈ 2.0^2 * 3.0 + 3.0^3
        @test grad[1] ≈ 2 * 2.0 * 3.0
        @test grad[2] ≈ 2.0^2 + 3 * 3.0^2
        @test M[1, 1] ≈ 2 * 3.0
        @test M[1, 2] ≈ 2 * 2.0
        @test M[2, 1] ≈ 2 * 2.0
        @test M[2, 2] ≈ 6 * 3.0
    end

    @testset "transcendental: sin(xy)" begin
        # ∇f = (y cos(xy), x cos(xy))
        # Hessian:
        #   ∂²/∂x² = -y² sin(xy)
        #   ∂²/∂x∂y = cos(xy) - xy sin(xy)
        #   ∂²/∂y² = -x² sin(xy)
        H = hessian(sin(x * y), [x, y])
        xv, yv = 1.0, 2.0
        _, _, M = H((xv, yv))
        xy = xv * yv
        @test M[1, 1] ≈ -yv^2 * sin(xy)
        @test M[1, 2] ≈ cos(xy) - xy * sin(xy)
        @test M[2, 1] ≈ cos(xy) - xy * sin(xy)
        @test M[2, 2] ≈ -xv^2 * sin(xy)
    end

    @testset "three variables: x² + y² + z² + xyz" begin
        H = hessian(x^2 + y^2 + z^2 + x*y*z, [x, y, z])
        _, _, M = H((1.0, 2.0, 3.0))
        # ∂²/∂x² = 2,  ∂²/∂y² = 2,  ∂²/∂z² = 2
        # ∂²/∂x∂y = z = 3,  ∂²/∂x∂z = y = 2,  ∂²/∂y∂z = x = 1
        @test M[1, 1] ≈ 2.0
        @test M[2, 2] ≈ 2.0
        @test M[3, 3] ≈ 2.0
        @test M[1, 2] ≈ 3.0; @test M[2, 1] ≈ 3.0
        @test M[1, 3] ≈ 2.0; @test M[3, 1] ≈ 2.0
        @test M[2, 3] ≈ 1.0; @test M[3, 2] ≈ 1.0
    end

    @testset "symmetry (Schwarz)" begin
        # ∂²f/∂x∂y = ∂²f/∂y∂x whenever f is C²
        for ex in (x^3 + 2*x^2*y + y^2, exp(x*y) + x^2*y, sin(x)*cos(y))
            H = hessian(ex, [x, y])
            _, _, M = H((0.7, -0.3))
            @test M[1, 2] ≈ M[2, 1]
        end
    end
end
