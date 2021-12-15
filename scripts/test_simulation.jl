using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using PyPlot
using Formatting

pygui(true)

##

function plotcrater(c::Crater)
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

t = 3.9
θₛ = π/4
rₑ = 1
Δ = 50
rmin = 100
nmax = 10_000_000
seed = 1

res = simulateimpacts(t, θₛ, rₑ, Δ, rmin=rmin, nmax=nmax, seed=seed)

##

#plot all the craters
for crater ∈ res.impactors
    plotcrater(crater)
end
#plot the extant shoreline segments
θ = [θₛ, θₛ]
for (ϕ₁,ϕ₂) ∈ res.segments
    plot([ϕ₁,ϕ₂], θ, color="C0")
end
#xlim(0, 2π)
ylim(0, π)
