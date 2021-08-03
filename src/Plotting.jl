using Diffusion
using StaticArrays
using Plots; gr()

function circleshape(x, y, r)
    θ = LinRange(0, 2*π, 500)
    x .+ r*sin.(θ), y .+ r*cos.(θ)
end

function plotreflectrow(rownumber, r)
    df = DataFrame(trial, :auto)
    rowx = []
    for x in 1:length(df[rownumber,:])
        push!(rowx, df[rownumber,:][x][1])
    end
    rowy = []
    for x in 1:length(df[rownumber,:])
        push!(rowy, df[rownumber,:][x][2])
    end
    scatter(rowx, rowy, legend = false)
    plot!(circleshape(0, 0, r), seriestype = [:shape], lw = 5, c = :red, linecolor = :lightgreen, fillalpha = 0.1, aspectratio = 2/2)
end

function animate(radius)
    anim = @animate for i = 1:length(trial[1])
        #plotreflectrow(i, radius)
        plotrow(i, radius)
    end 
    gif(anim, "graph.gif", fps = 5)
end

function animateforscript(trial, radius, savename)
    anim = @animate for i = 1:length(trial[1])
        #plotreflectrow(i, radius)
        plotrow(i, radius)
    end 
    gif(anim, savename, fps = 5)
end


function showrow(rownumber) 
    set = filter((x) -> length(x) >= rownumber, trial)
    setrow = []
    for z in 1:length(set)
        push!(setrow, set[z][rownumber])
    end
    return setrow
end

function plotrow(rownumber, r)
    setrow = showrow(rownumber)
    rowx = []
    for x in 1:length(setrow)
        push!(rowx, setrow[x][1])
    end
    rowy = []
    for x in 1:length(setrow)
        push!(rowy, setrow[x][2])
    end
    scatter(rowx, rowy, legend = false)
    plot!(circleshape(0, 0, r), seriestype = [:shape], lw = 5, c = :red, linecolor = :lightgreen, fillalpha = 0.1, aspectratio = 2/2)
end

trial = nmolreflect(10, 500, SA[10., 10.], 75, 1000)

animate(75)

