module Diffusion
using StaticArrays 

export absorb_boundary
export reflect_boundary 
export simmolreflect
export simmolabsorb
export nmolreflect
export nmolabsorb

###    boundary conditions     ###

function absorb_boundary(pos,r) #models an absorbing boundary-pos inside the boundary stay, but pos outside disappear 
    if(sum(abs2, pos) < r^2)
        return pos
    else
        return nothing
    end
end

function topolar(pos) #turns a cartesian coordinate into a polar coordinate 
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

function tocart(pos) #turns a polar coordinate into a cartesian coordinate
    x = pos[1]*cos(pos[2])
    y = pos[1]*sin(pos[2])
    return SA[x,y] 
end

function reflect(pos,r) #reflects a point radially over a circle 
    pos1 = topolar(pos)
    dist = pos1[1] - r
    pos2 = SA[r - dist,pos1[2]]
    posf = tocart(pos2)
    return posf
end

function reflect_boundary(pos,r) #models a reflective boundary. pos inside the boundary stay, pos outside the boundary are reflected back in
    if(sum(abs2, pos) <= r^2)
        return pos
    else
        pos = reflect(pos,r)
        return reflect_boundary(pos, r)
    end
end



###   simulating the movement of molecules inside a particular boundary   ### 

function movemol(pos, σ) #moves the molecule once  
    return pos + SA[σ*randn(),σ*randn()]
end

#= 
    simmol functions
    inputs: σ aka sqrt(2*d*deltat) aka the sd of the displacement, number of steps, initial position, radius
    outputs: vector of static arrays showing where the molecule has been
=#

function simmolreflect(σ, steps, pos, r) #simmol under refelctive boundary conditions 
    molhistory = [pos]
    for x in 1:steps
        pos = movemol(pos, σ)
        pos = reflect_boundary(pos, r)
        push!(molhistory, pos)
    end 
    return molhistory
end

function simmolabsorb(σ, steps, pos, r) #simmol under absorbing boundary conditions
    molhistory = [pos]
    for x in 1:steps
        pos = movemol(pos, σ)
        pos = absorb_boundary(pos, r)
        if pos === nothing
            return molhistory
        else
        push!(molhistory, pos)
        end
    end 
    return molhistory
end

function nmolreflect(σ, steps, pos, r, molecules) #simulates multiple molecules at once-reflective
    molhistories = []
    for x in 1:molecules
        push!(molhistories, simmolreflect(σ, steps,pos,r))
    end
    return molhistories
end

function nmolabsorb(σ, steps, pos, r, molecules) #simulates multiple molecules at once-absorbing
    molhistories = []
    for x in 1:molecules
        push!(molhistories, simmolabsorb(σ, steps,pos,r))
    end
    return molhistories
end

end
