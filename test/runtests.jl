using Diffusion
using Test

@testset "Diffusion.jl" begin
    # Write your tests here.
    @test absorb_boundry(1,2) == 3
end

