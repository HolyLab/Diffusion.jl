using Diffusion
using Test
using StaticArrays

@testset "Diffusion.jl" begin
    # Write your tests here.
    @test absorb_boundry(SA[3,4],6) == SA[3,4]
    @test absorb_boundry(SA[3,4],5) === nothing
    @test absorb_boundry(SA[3,4],4) === nothing 

    @test reflect_boundary(SA[3,4],6) == SA[3,4]
    
end

