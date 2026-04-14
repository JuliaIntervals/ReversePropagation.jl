using IntervalArithmetic, IntervalArithmetic.Symbols
using IntervalBoxes

eq(a, b) = isequal_interval(bareinterval(a), bareinterval(b))
eq(a::Tuple, b::Tuple) = all(eq.(a, b))

@testset "forward_backward_contractor" begin

    vars = @variables x, y
    @variables a  # parameter

    C = forward_backward_contractor(x + 1, x)
    @test eq(C(IntervalBox(-10..10), 2..3)[1], (1..2, ))


    ex = x^2 + y^2

    C = forward_backward_contractor(ex, vars)
    @test eq(C(IntervalBox(-10..10, -10..10), 0..1), ( (-1..1, -1..1), 0..200 ))
end

@testset "forward_backward_contractor with parameter" begin

    vars = @variables x, y
    @variables a  # parameter

    ex = x + a

    C = forward_backward_contractor(ex, [x], [a])

    @test eq(C(IntervalBox(-10..10), 5..6, 2..3), ( (2..4, ), (-8..13) ) )


    ex = x^2 + a * y^2

    C = forward_backward_contractor(ex, vars, [a])
    @test eq(C(IntervalBox(-10..10, -10..10), 0..1, 1..1), ( (-1..1, -1..1), 0..200 ))
end

@testset "SSAFunction with forward_backward_contractor" begin

    vars = @variables x, y

    ex = x^2 + y^2

    # Build SSAFunction manually via cse_equations, then pass to forward_backward_contractor
    ssa = ReversePropagation.cse_equations(ex)
    @test ssa isa SSAFunction
    @test length(ssa.code) > 0

    C = forward_backward_contractor(ssa, vars)
    @test eq(C(IntervalBox(-10..10, -10..10), 0..1), ( (-1..1, -1..1), 0..200 ))
end

@testset "derivative_contractor" begin

    vars = @variables x, y

    @testset "simple polynomial" begin
        # d/dx(x^2 + y^2) = 2x
        ex = x^2 + y^2
        C = derivative_contractor(ex, x, vars)

        # 2x == 0 → x = 0
        r = C(IntervalBox(-10..10, 2), interval(0, 0))
        @test eq(r[1], (0..0, -10..10))

        # 2x ∈ [-2, 2] → x ∈ [-1, 1]
        r2 = C(IntervalBox(-10..10, 2), interval(-2, 2))
        @test eq(r2[1], (-1..1, -10..10))
    end

    @testset "with transcendental functions" begin
        # d/dx(sin(x*y) + x^2) = y*cos(x*y) + 2x
        ex = sin(x*y) + x^2
        C = derivative_contractor(ex, x, vars)

        # On [-2,2]^2: contractor should run without error and contract
        r = C(IntervalBox(-2..2, 2), interval(0, 0))
        # x should be contracted (derivative = 0 is a strict constraint)
        @test all(x -> diam(x) <= 4, r[1])
    end

    @testset "derivative w.r.t. second variable" begin
        # d/dy(x^2 + y^2) = 2y
        ex = x^2 + y^2
        C = derivative_contractor(ex, y, vars)

        # 2y == 0 → y = 0
        r = C(IntervalBox(-10..10, 2), interval(0, 0))
        @test eq(r[1], (-10..10, 0..0))
    end

    @testset "derivative_ssa" begin
        ex = x^2 + y^2
        ssa = derivative_ssa(ex, x, vars)
        @test ssa isa SSAFunction
        @test length(ssa.code) > 0
    end
end

@testset "bare intervals" begin

    vars = @variables x, y
    @variables a  # parameter

    ex = x^exact(2) + y^exact(2)

    C = forward_backward_contractor(ex, vars)
    X = IntervalBox(bareinterval(-10..10), 2)
    constraint = bareinterval(0..1)

    @test eq(C(X, constraint), ( (bareinterval(-1..1), bareinterval(-1..1)), bareinterval(0..200) ) )


    # with parameter:
    # ex = x^exact(2) + a * y^exact(2)

    # C = forward_backward_contractor(ex, vars, [a])

    # X = IntervalBox(bareinterval(-10..10), 2)
    # constraint = bareinterval(0..1)

    # @test eq(C(X, constraint, [bareinterval(1..1)]),
    #     ( (bareinterval(-1..1), bareinterval(-1..1)), bareinterval(0..200) ) )

end