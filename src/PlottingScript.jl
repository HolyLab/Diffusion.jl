using Diffusion
using Unitful: cm, s
using StaticArrays

#snapshots
#how many snapshots do you want to take?
numberofsnapshots = 3
#what time steps do you want snapshots of?
snapshottimesteps = [250, 500, 1000]

#trial name
name = "Trial 1 Amyl Acetate"
#length of each time step
t = 0.02 
#diffusion constant in air
D = 0.0067
#number of time steps
steps = 999
#original position
pos = SA[0.0,0.0]
#radius of the circle
r = 0.4
#number of simulated molecules
mols = 10000
#fast sniff time
fastsniff = 250
#middle sniff time
middlesniff = 500
#regular respiration time 
respiration = 1000
#Standard Deviation of the Displacement
σ = sqrt(2*D*t)

#create simulation data
trialreflect = nmolreflect(σ, steps, pos, r, mols)
trialabsorb = nmolabsorb(σ, steps, pos, r, mols)

#create analysis 
reflectavg = molavgbytime(trialreflect)
reflectorigindist = []
for i in 1:length(reflectavg)
    push!(reflectorigindist, sqrt(sum(abs2, reflectavg[i])))
end
absorbavg = molavgbytime(trialabsorb)
absorborigindist = []
for i in 1:length(absorbavg)
    push!(absorborigindist, sqrt(sum(abs2, absorbavg[i])))
end
reflectvarcoor = molvarbytime(trialreflect)
reflectvar = []
for i in 1:length(reflectvarcoor)
    push!(reflectvar, sqrt(sum(abs2, reflectvarcoor[i])))
end
absorbvarcoor = molvarbytime(trialabsorb)
absorbvar = []
for i in 1:length(absorbvarcoor)
    push!(absorbvar, sqrt(sum(abs2, absorbvarcoor[i])))
end
reflectnumber = numberremaining(molpathlengths(trialreflect))
absorbnumber = numberremaining(molpathlengths(trialabsorb))
reflectpercent = percentremaining(molpathlengths(trialreflect))
absorbpercent = percentremaining(molpathlengths(trialabsorb))

#create visualizations
#animations
trial = trialreflect
animateforscript(trialreflect, r, "$name reflect animation.gif")
trial = trialabsorb
animateforscript(trialabsorb, r, "$name absorb animation.gif")

#plots
percentplot = plot(1:length(reflectpercent), reflectpercent, color = :blue, xlabel = "time", ylabel = "% remaining", title = "$name % of Mols Remaining vs Time", label = "reflect", linewidth = 5, grid = false, legend = :right)
plot!(1:length(absorbpercent), absorbpercent, color = :red, label = "absorb", linewidth = 5)
vline!([fastsniff], label = "fast sniff", color = :magenta, linewidth = 5)
vline!([middlesniff], label = "middle sniff", color = :lightgreen, linewidth = 5)
vline!([respiration], label = "regular respiration", color = :black, linewidth = 5)
png(percentplot, "$name percent remaining.png")

varianceplot = plot(1:length(reflectvar), reflectvar, color = :blue, xlabel = "time", ylabel = "variance", title = "$name Variance vs Time", label = "reflect", linewidth = 5, grid = false, legend = :right)
plot!(1:length(absorbvar), absorbvar, color = :red, label = "absorb", linewidth = 5)
vline!([fastsniff], label = "fast sniff", color = :magenta, linewidth = 5)
vline!([middlesniff], label = "middle sniff", color = :lightgreen, linewidth = 5)
vline!([respiration], label = "regular respiration", color = :black, linewidth = 5)
png(varianceplot, "$name variance.png")

distplot = plot(1:length(reflectorigindist), reflectorigindist, color = :blue, xlabel = "time", ylabel = "distance to origin", title = "$name dist to origin vs time", label = "reflect", linewidth = 5, grid = false, legend = :right)
plot!(1:length(absorborigindist), absorborigindist, color = :red, label = "absorb", linewidth = 5)
vline!([fastsniff], label = "fast sniff", color = :magenta, linewidth = 5)
vline!([middlesniff], label = "middle sniff", color = :lightgreen, linewidth = 5)
vline!([respiration], label = "regular respiration", color = :black, linewidth = 5)
png(distplot, "$name dist to origin.png")

#snapshots
trial = trialreflect
for i in 1:numberofsnapshots
    rowname = snapshottimesteps[i]
    snapplot = plotrow(snapshottimesteps[i], r)
    scatter!([pos[1]], [pos[2]], legend = false, markershape = :cross, markersize = 7)
    png(snapplot, "$name Snapshot of Timestep $rowname Reflect.png")
end
trial = trialabsorb
for i in 1:numberofsnapshots
    rowname = snapshottimesteps[i]
    snapplot = plotrow(snapshottimesteps[i], r)
    scatter!([pos[1]], [pos[2]], legend = false, markershape = :cross, markersize = 7)
    png(snapplot, "$name Snapshot of Timestep $rowname Absorb.png")
end
