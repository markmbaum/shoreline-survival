using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using UnPack
using ShorelineSurvival

include(scriptsdir("plotting_functions.jl"))

##

t = 4.1
θₛ = 2π/5
rmin = 1000
nmax = 1000
seed = 5
R = ♂ᵣ

##

figure()
#plot the shoreline ring
θ = fill(θₛ, 100)
ϕ = LinRange(0, π, 100)
x, y, z = sph2cart(θ, ϕ, R)
plot3D(x, y, z, color="C3", linewidth=0.75)
#plot all the craters
res = simulateimpacts(t, θₛ, rmin=rmin, nmax=nmax, seed=seed)
for crater ∈ GlobalPopulation(t, rmin=rmin, nmax=nmax, seed=seed)
    @unpack θ, ϕ, r = crater
    if 0.05 < ϕ < π - 0.05
        if abs(θ - θₛ) < r/R
            plotcirc(θ, ϕ, r, R=R, color=:k, linewidth=0.75)
        else
            plotcirc(θ, ϕ, r, R=R, color=:gray, linewidth=0.5)
        end
    end
end
θ = fill(θₛ, 25)
for seg ∈ res.segments
    if (0.05 < seg[1] < π - 0.05) & (0.05 < seg[2] < π - 0.05)
        ϕ = LinRange(seg[1], seg[2], 25)
        x, y, z = sph2cart(θ, ϕ, R)
        plot3D(x, y, z, color="C0", linewidth=0.75)
    end
end
axis("off")
gca()[:view_init](0,90)
tight_layout()

savefig(plotsdir("example_simulation_isolatitude.png"), dpi=1200)

##

segments = readsegments(
    datadir(
        "exp_pro",
        "parker_1989_contact_1a.csv"
    ),
    minarc=0.01
);

figure()
#plot the shoreline segments
for s ∈ segments
    if (π <= s.a.ϕ <= 2π) & (π <= s.b.ϕ <= 2π)
        plotseg(s, R=R)
    end
end
#plot all the craters
P = GlobalPopulation(t, rmin=rmin, nmax=nmax, seed=seed)
res = simulateimpacts(t, segments, rmin=rmin, nmax=nmax, seed=seed)
foreach(plotseg, segments)
for s ∈ segments
    if (π <= s.a.ϕ <= 2π) & (π <= s.b.ϕ <= 2π)
        plotseg(s, R=R, color=:C3)
    end
end
for s ∈ res.segments
    if (π <= s.a.ϕ <= 2π) & (π <= s.b.ϕ <= 2π)
        plotseg(s, R=R, color=:C0)
    end
end
for crater ∈ P
    #take only one hemisphere
    @unpack θ, ϕ, r = crater
    if π + 0.05 < ϕ < 2π - 0.05
        if crater ∈ res.impactors
            plotcirc(θ, ϕ, r, R=R, color=:k, linewidth=0.75)
        else
            plotcirc(θ, ϕ, r, R=R, color=:gray, linewidth=0.5)
        end
    end
end
axis("off")
gca()[:view_init](0,90)
tight_layout()

savefig(plotsdir("example_simulation_mapped.png"), dpi=1200)
