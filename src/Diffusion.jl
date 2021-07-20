module Diffusion
using StaticArrays 

export absorb_boundary
export reflect_boundary 

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
    if(sum(abs2, pos) < r^2)
        return pos
    else
        return reflect(pos,r)
    end
end

end
