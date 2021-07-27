using Diffusion
using Test
using StaticArrays

@testset "Diffusion.jl" begin
    # tests for absorbing and reflective boundaries
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
    
    #tests for simulating the movement of molecules
    @test length(simmolreflect(1, 4, SA[0.,0.], 3)) == 5
    @test typeof(simmolreflect(1, 4, SA[1.,1.], 4)) == Vector{SVector{2, Float64}}
    @test length(simmolabsorb(1, 4, SA[0.,0.], 3)) <= 5
    @test typeof(simmolabsorb(1, 4, SA[1.,1.], 4)) == Vector{SVector{2, Float64}}
    @test length(nmolreflect(1, 4, SA[1., 1.], 3, 100)) == 100
    @test typeof(nmolreflect(1, 4, SA[1., 1.], 3, 100)) == Vector{Any}
    @test length(nmolabsorb(1, 4, SA[1., 1.], 3, 100)) == 100
    @test typeof(nmolabsorb(1, 4, SA[1., 1.], 3, 100)) == Vector{Any}
end

