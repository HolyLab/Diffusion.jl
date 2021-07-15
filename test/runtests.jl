using Diffusion
using Test
using StaticArrays

@testset "Diffusion.jl" begin
    # Write your tests here.
    @test absorb_boundry(SA[3,4],6) == SA[3,4]
    @test absorb_boundry(SA[3,4],5) === nothing
    @test absorb_boundry(SA[3,4],4) === nothing 
    @test slope(SA[1,2],SA[3,4]) == 1
    @test midpoint(SA[1,1],SA[3,3]) == SA[2,2]
    @test isapprox(circleintersect(SA[1,0], SA[7,0], 2),SA[2,0])
    @test rslope(SA[2,3]) == 3/2
    @test perpslope(3) == -1/3
    @test solve(1,-2,0,3) == SA[1,1]
    @test posFinal(SA[1,2],SA[1.5,1.5]) == SA[2,1]
    @test isapprox(reflect_boundary(SA[0,2], SA[4,2], sqrt(8)), SA[2,0])
    
end

