using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using PyPlot
using Formatting
using ProfileView
using BenchmarkTools

pygui(true)

##

function plotcrater(c::GlobalCrater)
    x, y, z = sphcirc(c, N=50)
    θ, ϕ, _ = cart2sph(x, y, z)
    if any(ϕ .> 7π/4) & any(ϕ .< π/4)
        b = ϕ .< π
        plot(ϕ[b], θ[b], color="k")
        b = ϕ .> π
        plot(ϕ[b], θ[b], color="k")
    else
        plot(ϕ, θ, color="k")
    end
    nothing
end

##

t = 4
θₛ = π/4
rₑ = 1.5
Δ = 50
rmin = 100
nmax = 1e8
seed = 1

#ProfileView.@profview begin
    res = globalsimulation(
        t,
        θₛ,
        rₑ,
        Δ,
        rmin=rmin,
        nmax=nmax,
        seed=seed
    );
    print(res);
#end;

##

@btime begin
    globalsimulation(
        4,
        π/4,
        1,
        0,
        rmin=100,
        nmax=1e7,
        seed=1,
        show=false
    );
end;

##

if nmax > 1e4
    error("Too many craters!")
else
    figure()
    #plot all the craters
    for crater ∈ res.impactors
        plotcrater(crater)
    end
    #plot the extant shoreline segments
    θ = [θₛ, θₛ]
    for (ϕ₁,ϕ₂) ∈ res.segments
        plot([ϕ₁,ϕ₂], θ, color="C0")
    end
    gca()[:axis]("square")
end