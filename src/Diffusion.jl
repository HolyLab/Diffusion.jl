module Diffusion
using StaticArrays 

export absorb_boundary
export reflect_boundary 
export simmolreflect
export simmolabsorb

function absorb_boundary(pos,r)
    if(sum(abs2, pos) < r^2)
        return pos
    else
        return nothing
    end
end

function topolar(pos)
    if pos[1] >= 0 && pos[2] >= 0
        theta = atan(pos[2]/pos[1])
    elseif pos[1] >=0 && pos[2] <= 0
        theta = atan(pos[2]/pos[1]) + 2π
    else
        theta = atan(pos[2]/pos[1]) + π
    end
    r = sqrt(sum(abs2, pos))
    return SA[r,theta]
end

function tocart(pos)
    x = pos[1]*cos(pos[2])
    y = pos[1]*sin(pos[2])
    return SA[x,y] 
end

function reflect(pos,r)
    pos1 = topolar(pos)
    dist = pos1[1] - r
    pos2 = SA[r - dist,pos1[2]]
    posf = tocart(pos2)
    return posf
end

function reflect_boundary(pos,r)
    if(sum(abs2, pos) <= r^2)
        return pos
    else
        pos = reflect(pos,r)
        return reflect_boundary(pos, r)
    end
end

#= create function simmol
    inputs: movement per step, number of steps, initial position, radius
    outputs: vector of static arrays (showing where the molecule has been)
=#

function movemol(pos) #does not actually model movement of a particle yet
    return pos + SA[rand(-3:3),rand(-3:3)]
end

function simmolreflect(steps,pos,r) #currently does not model movement per step correctly
    molhistory = [pos]
    for x in 1:steps
        pos = movemol(pos)
        pos = reflect_boundary(pos, r)
        push!(molhistory, pos)
    end 
    return molhistory
end

function simmolabsorb(steps,pos,r) #currently does not model movement per step correctly
    molhistory = [pos]
    for x in 1:steps
        pos = movemol(pos)
        pos = absorb_boundary(pos, r)
        if pos === nothing
            return molhistory
        else
        push!(molhistory, pos)
        end
    end 
    return molhistory
end

end
