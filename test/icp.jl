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