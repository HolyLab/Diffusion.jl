module Diffusion
using StaticArrays 
using CoordinateTransformations 

export absorb_boundary
export reflect_boundary 
export simmolreflect
export simmolabsorb
export nmolreflect
export nmolabsorb
export molpathlengths
export numberremaining
export molavgbytime
export molvarbytime

###    boundary conditions     ###

function absorb_boundary(pos,r) #models an absorbing boundary-pos inside the boundary stay, but pos outside disappear 
    if(sum(abs2, pos) < r^2)
        return pos
    else
        return nothing
    end
end

function topolar(pos) #turns a cartesian coordinate into a polar coordinate 
    #if pos[1] >= 0 && pos[2] >= 0
        theta = atan(pos[2], pos[1])
    #elseif pos[1] >=0 && pos[2] <= 0
    #    theta = atan(pos[2]/pos[1]) + 2π
    #else
    #    theta = atan(pos[2]/pos[1]) + π
    #end
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
    molhistories = Vector{typeof(pos)}[]
    for x in 1:molecules
        push!(molhistories, simmolreflect(σ, steps,pos,r))
    end
    return molhistories
end

function nmolabsorb(σ, steps, pos, r, molecules) #simulates multiple molecules at once-absorbing
    molhistories = Vector{typeof(pos)}[]
    for x in 1:molecules
        push!(molhistories, simmolabsorb(σ, steps,pos,r))
    end
    return molhistories
end

### Analysis ###
 
function molremaining(σ, steps, pos, r, molecules, boundary)
    if boundary == "absorb"
        trial = nmolabsorb(σ, steps, pos, r, molecules)
    elseif boundary == "reflect"
        trial = nmolreflect(σ, steps, pos, r, molecules)
    else
        throw(DomainError("Please specify boudary condition as 'absorb' or 'reflect'"))
    end
    lengths = molpathlengths(trial)
    return percentremaining(lengths)
    end
    
    function molpathlengths(trial) #input the molhistories given my nmolreflect or nmolabsorb
    lengthsofmolpaths = []         #outputs a vector showing the lengths of each molecules' path
    for x in 1:length(trial)
        push!(lengthsofmolpaths, length(trial[x]))
    end
    return lengthsofmolpaths
    end
    
    
    function numberremaining(lengthsofmolpaths) #input a vector containing the length of each molecules path
        numberremain = []                       #out puts a vector containing the number of molecules remaining each time step
            for x in 1:maximum(lengthsofmolpaths)
                number = 0.0
                for y in 1:length(lengthsofmolpaths)
                    if lengthsofmolpaths[y] >= x 
                        number = number + 1
                    end
                end
                push!(numberremain, number)
            end
        return numberremain
    end
    
    function percentremaining(lengthsofmolpaths)                                   #input a vector containing the length of each molecules path
    percentremain = numberremaining(lengthsofmolpaths)/length(lengthsofmolpaths)   #out puts a vector containing the percent of molecules remaining each time step 
    return percentremain
    end
    
    #average and variance of position of all remaining molecules as a function of time. Input simulation specification, output vector of means or variances by time step
    
    function molavgbytime(σ, steps, pos, r, molecules, boundary) 
        if boundary == "absorb"
            trial = nmolabsorb(σ, steps, pos, r, molecules)
        elseif boundary == "reflect"
            trial = nmolreflect(σ, steps, pos, r, molecules)
        else
            throw(DomainError("Please specify boudary condition as 'absorb' or 'reflect'"))
        end
        lengths = molpathlengths(trial)
        means = []
        for row in 1:maximum(lengths)
            set = filter((x) -> length(x) >= row, trial)
            setrow = []
            for z in 1:length(set)
                push!(setrow, set[z][row])
            end
            avg = mean(setrow)
            push!(means, avg)
        end
        return means
    end
    
    
    function molvarbytime(σ, steps, pos, r, molecules, boundary) 
        if boundary == "absorb"
            trial = nmolabsorb(σ, steps, pos, r, molecules)
        elseif boundary == "reflect"
            trial = nmolreflect(σ, steps, pos, r, molecules)
        else
            throw(DomainError("Please specify boudary condition as 'absorb' or 'reflect'"))
        end
        lengths = molpathlengths(trial)
        vars = []
        for row in 1:maximum(lengths)
            set = filter((x) -> length(x) >= row, trial)
            setrow = []
            for z in 1:length(set)
                push!(setrow, set[z][row])
            end
            variance = var(setrow)
            push!(vars, variance)
        end
        return vars
    end
    

end
