using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using Plots
using BenchmarkTools

##

function testintersection(θₛ, r)
    #get a crater that crosses θₛ
    c = Crater(r)
    while abs(c.θ - θₛ)*♂ᵣ > r
        c = Crater(r)
    end
    #find intersection interval
    ϕ₁, ϕ₂ = intersection(c, θₛ)
    #draw an intersection line
    θ = fill(θₛ, 10)
    if ϕ₁ > ϕ₂
        ϕ = [LinRange(0, ϕ₁, 5); LinRange(ϕ₂, 2π, 5)]
    else
        ϕ = LinRange(ϕ₁, ϕ₂, 10)
    end
    #draw the crater
    x, y, z = sphcirc(c)
    Θ, Φ, _ = cart2sph(x, y, z)
    #return for plotting by whatever means
    return θ, ϕ, Θ, Φ
end

##

θₛ = 2π/4
r = 1e5
θ, ϕ, Θ, Φ = testintersection(θₛ, r)
p = plot(; legend=false)
plot!(p, ϕ, θ, markersize=2, color=:orange)
plot!(p, Φ, Θ, markersize=2, color=:red)

##

function timeintersection(θₛ, r)
    #get a crater that crosses θₛ
    c = Crater(r)
    while abs(c.θ - θₛ)*♂ᵣ > r
        c = Crater(r)
    end
    #time the intersection
    @btime intersection($c, $θₛ)
end