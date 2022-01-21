using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using MultiAssign
using PyPlot
using ProfileView
using BenchmarkTools

pygui(true)

## handy 2D plotting functions

function plotcrater(c::GlobalCrater, color="k", linewidth=1; N=50)
    x, y, z = sphcirc(c, N=N)
    θ, ϕ, _ = cart2sph(x, y, z)
    if any(ϕ .> 7π/4) & any(ϕ .< π/4)
        b = ϕ .< π
        plot(ϕ[b], θ[b], color=color, linewidth=linewidth)
        b = ϕ .> π
        plot(ϕ[b], θ[b], color=color, linewidth=linewidth)
    else
        plot(ϕ, θ, color=color, linewidth=linewidth)
    end
    nothing
end

function plotsegment(s::SphericalSegment)
    plot(
        [s.a.ϕ, s.b.ϕ],
        [s.a.θ, s.b.θ],
        color='C'*string(rand(1:8)),
        alpha=0.9,
        linewidth=1,
        zorder=10
    )
    nothing
end

function plotgreatcircle(C::GreatCircle)
    t = LinRange(0, 𝛕, 1000)
    @multiassign θ, ϕ = zeros(length(t))
    for i = 1:length(t)
        g = C(t[i])
        θ[i], ϕ[i] = cart2usph(g...)
    end
    plot(ϕ, θ, "k.", markersize=1, alpha=0.5)
end

function plotgreatcircle(s::SphericalSegment)
    plotgreatcircle(GreatCircle(s))
end

## mapped putative shoreline coordinates

segments = readsegments(
    datadir("exp_pro", "parker_1989_contact_1a.csv"),
    minarc=0.04
);

## for testing with a straight line around the equator

ϕ = LinRange(0, 2π, 50)
segments = [SphericalSegment((π/2, ϕ[i]), (π/2, ϕ[i+1])) for i ∈ 1:length(ϕ)-1];

##

t = 4
rₑ = 1
Δ = 0
rmin = 1e2
nmax = 1e6
seed = 1

#ProfileView.@profview begin
    res = globalsimulation(
        t,
        segments,
        rₑ,
        Δ,
        rmin=rmin,
        nmax=nmax,
        seed=seed,
        show=true
    )
#end;

println("smallest gap between segments: ", minimum(gapdistances(res)))

##

@btime begin
    globalsimulation(
        4,
        $segments,
        1,
        0,
        rmin=100,
        nmax=1e6,
        seed=1
    )
end;

##

figure()
impactors = Set(res.impactors)
for crater ∈ GlobalPopulation(t, rmin=max(rmin,Δ), nmax=nmax, seed=seed)
    crater *= rₑ
    if crater ∈ impactors
        plotcrater(crater, "r", 1, N=500)
    else
        plotcrater(crater, "k", 1, N=50)
    end
end
foreach(plotsegment, res.segments)
#plotgreatcircle.(segments)
#xlim(0, 𝛕)
#ylim(0, π)
