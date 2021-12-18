using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using LinearAlgebra: ⋅
using MultiAssign
using PyPlot
using ProfileView
using BenchmarkTools

pygui(true)

##

function plotcrater(c::Crater, color="k", linewidth=1)
    x, y, z = sphcirc(c, N=50)
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

function plotsegments(segments::Vector{SphericalSegment})
    for s ∈ segments
        plot(
            [s.a.ϕ, s.b.ϕ],
            [s.a.θ, s.b.θ],
            color="k",
            alpha=0.9,
            linewidth=2,
            zorder=10
        )
    end
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

##

#shoreline coordinates
fn = datadir("exp_pro", "parker_1989_contact_1a.csv")
#read the coordinates into segments with appropriate spacing
segments = readsegments(fn, minarc=0.02)#0.02)
#segments = rand(segments, 5)

##

t = 4
θₛ = π/4
rₑ = 1
Δ = 0
rmin = 100
nmax = 1e3
seed = 1

#ProfileView.@profview begin
    res = simulateimpacts(
        t,
        segments,
        rₑ,
        Δ,
        rmin=rmin,
        nmax=nmax,
        seed=seed,
        show=true
    )
   println(res)
#end;

##

@btime begin
    simulateimpacts(
        $t,
        $segments,
        $rₑ,
        $Δ,
        rmin=100,
        nmax=1e3,
        seed=$seed,
        show=false
    )
end

##

figure()
for crater ∈ GlobalPopulation(t, rmin=max(rmin,Δ), nmax=nmax, seed=seed)
    crater *= rₑ
    if crater ∈ res.impactors
        plotcrater(crater, "r", 1)
    else
        plotcrater(crater, "k", 1)
    end
end
plotsegments(res.segments)
#plotgreatcircle.(segments)
xlim(0, 𝛕)
ylim(0, π)