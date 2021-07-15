module Diffusion
using StaticArrays 

export absorb_boundry
export reflect_boundary 
export slope
export midpoint
export circleintersect
export rslope
export perpslope
export solve
export posFinal

function absorb_boundry(pos,r)
    if(((pos[1]^2)+(pos[2]^2))^(1/2) < r)
    return pos
    else
    return nothing
    end
end

function slope(pos1,pos2)               #finds the slope between two points
    return (pos1[2]-pos2[2])/(pos1[1]-pos2[1])
end

function midpoint(pos1,pos2)            #finds the midpoint between two points
    x = (pos1[1]+pos2[1])/2
    y = (pos1[2]+pos2[2])/2
    return SA[x,y]
 end

function circleintersect(pos1,pos2,r)   #estimates the intersection between the line formed by two points and a circle of radius r
    mp = midpoint(pos1,pos2)
    if(isapprox(((pos1[1]^2)+(pos1[2]^2))^(1/2),r))
        return pos1
    elseif(((mp[1]^2)+(mp[2]^2))^(1/2) < r) 
        circleintersect(mp, pos2, r)
    else 
        circleintersect(pos1, mp, r)
    end
end

function rslope(posC)
    return posC[2]/posC[1]
end

function perpslope(m)
    return -1/m
end

function solve(m1,m2,x1,y1) #slope of radius, slope of perpendicular, point on perpendicular line (x and y)
    A = [-m1 1;-m2 1]
    b = [0, y1 - m2*x1]
    temp = A\b
    return SA[temp[1],temp[2]]
end

function posFinal(pos1,posM)
    deltax = posM[1] - pos1[1]
    deltay = posM[2] - pos1[2]
    return SA[posM[1] + deltax, posM[2] + deltay]
end

function reflect_boundary(pos1,pos2,r)
    if(((pos2[1]^2)+(pos2[2]^2))^(1/2) < r)
    return pos
    else
        posC = circleintersect(pos1,pos2,r)
        mr = rslope(posC)
        mp = perpslope(mr)
        posM = solve(mr,mp,pos1[1],pos1[2])
        return posFinal(pos1,posM)
    end
end

end
