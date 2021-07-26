using Diffusion
using Test
using StaticArrays

@testset "Diffusion.jl" begin
    # Write your tests here.
    @test absorb_boundary(SA[3,4],6) == SA[3,4]
    @test absorb_boundary(SA[3,4],5) === nothing
    @test absorb_boundary(SA[3,4],4) === nothing 
    @test reflect_boundary(SA[3,4],6) == SA[3,4]
    @test reflect_boundary(SA[3,4],5) ≈ SA[3,4]
    @test reflect_boundary(SA[4,4],sqrt(18)) ≈ SA[2,2] 
    
    for r in 0.5:0.1:1.2, θ in 0:π/8:2π
        pt = SA[r*cos(θ), r*sin(θ)]
        rbound = 0.8
        rreflect = 2*rbound - r
        @test reflect_boundary(pt, rbound) ≈ (r <= rbound ? pt : SA[rreflect*cos(θ), rreflect*sin(θ)])
    end
    
    @test length(simmolreflect(4, SA[0.,0.], 3)) == 5
    @test typeof(simmolreflect(4, SA[1.,1.], 4)) == Vector{SVector{2, Float64}}
    @test length(simmolabsorb(4, SA[0.,0.], 3)) <= 5
    @test typeof(simmolabsorb(4, SA[1.,1.], 4)) == Vector{SVector{2, Float64}}
end

