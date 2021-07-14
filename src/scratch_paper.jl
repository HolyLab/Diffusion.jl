using StaticArrays

function slope(pos1,pos2)               #finds the slope between two points
    return (pos1[2]-pos2[2])/(pos1[1]-pos2[1])
end

function midpoint(pos1,pos2)            #finds the midpoint between two points
    x = (pos1[1]+pos2[1])/2
    y = (pos1[2]+pos2[2])/2
    return SA[x,y]
 end

function distance(pos1,pos2)            #finds the distance between two points
    return ((pos2[1]-pos1[1])^2+(pos2[2]-pos1[2])^2)^(1/2)
end

function circleintersect(pos1,pos2,r)   #estimates the intersection between the line formed by two points and a circle of radius r
    if(!(((pos1[1]^2)+(pos1[2]^2))^(1/2) < r) || !(((pos2[1]^2)+(pos2[2]^2))^(1/2) > r))
throw(DomainError("pos1 must be inside the circle and pos2 must be outside the circle"))
    end
    mp = midpoint(pos1,pos2)
    if(isapprox(((pos1[1]^2)+(pos1[2]^2))^(1/2),r))
        return pos1
    elseif(((mp[1]^2)+(mp[2]^2))^(1/2) < r) 
        circleintersect(mp, pos2, r)
    else 
        circleintersect(pos1, mp, r)
    end
end

function circleslope(pos1)              #finds the slope of a circle at a point
    return -pos1[1]/pos1[2]
end

function angleofincidence(pos1,pos2,r)  #finds the angle of incidence of the line between pos1 and pos 2 as it reflects of a circle with radius r
    posC = circleintersect(pos1,pos2,r)
    m1 = slope(pos1,pos2)
    m2 = circleslope(posC)
    return atan((m1-m2)/(1+m1*m2))
end

function returnslope(pos1,pos2,r)       #finds the slope of the line between pos1 and pos 2 after reflection on a circle with radius r
    theta = angleofincidence(pos1,pos2,r) 
    posC = circleintersect(pos1, pos2, r)
    m2 = circleslope(posC)
    return (tan(-theta)+m2)/(1-tan(-theta)*m2)
end

function finallocation(pos1,pos2,r)     #currently broken-tries to find where a particle moving from pos1 to pos2 would reflect to after hitting a circle with radius r
    m = returnslope(pos1,pos2,r)
    posC = circleintersect(pos1,pos2,r)
    dist = distance(posC,pos2)
    posTemp1 = SA[posC[1] + 100000, posC[2] + 100000*m]
    posFinalA = circleintersect(posC, posTemp1, dist)       #the broken part is here-circle intersect uses a circle centered on 0, I need a circle centered on posC
    posTemp2 = SA[posC[1] - 100000, posC[2] - 100000*m]
    posFinalB = circleintersect(posC, posTemp2, dist)
    if(((posFinalA[1]^2)+(posFinalA[2]^2))^(1/2) < r)
        return posFinalA
    else
        return posFinalB
    end
end 

using Test

@testset "Diffusion.jl" begin
    @test slope(SA[1,2],SA[3,4]) == 1
    @test midpoint(SA[1,1],SA[3,3]) == SA[2,2]
    @test distance(SA[1,1],SA[2,1]) == 1
    @test isapprox(circleintersect(SA[1,0], SA[7,0], 2),SA[2,0])
    @test circleslope(SA[1,2]) == -0.5 
    @test angleofincidence(SA[1,1], SA[3,3], 2) == 90π/180
    @test isapprox(returnslope(SA[1,1], SA[3,3], 2), 1)
    #I don't know how to test final location
end
