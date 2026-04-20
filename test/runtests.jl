using ReversePropagation
using Symbolics
using Test

@testset "ReversePropagation.jl" begin

    include("gradient.jl")
    include("icp.jl")
    include("hardening.jl")
    include("tangent.jl")

end

