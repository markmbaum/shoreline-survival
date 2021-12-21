using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using BenchmarkTools

include(scriptsdir("plotting_functions.jl"))

##

function intersection(θ, ϕ, r, Θ, R)
    C = cos(r/R)
    S = sin(Θ)
    z = cos(Θ)
    a, b, c = sph2cart(θ, ϕ)
    
    x₁ = (a*C - a*c*z - sqrt(-b^2*C^2 + a^2*b^2*S^2 + b^4*S^2 + 2*b^2*c*C*z - b^2*c^2*z^2))/(a^2 + b^2)

    y₁ = (C - (a^2*C)/(a^2 + b^2) - c*z + (a^2*c*z)/(a^2 + b^2) + (a*sqrt(-b^2*(C^2 - a^2*S^2 - b^2*S^2 - 2*c*C*z + c^2*z^2)))/(a^2 + b^2))/b

    x₂ = (a*C - a*c*z + sqrt(-b^2*C^2 + a^2*b^2*S^2 + b^4*S^2 + 2*b^2*c*C*z - b^2*c^2*z^2))/(a^2 + b^2)

    y₂ = (C - (a^2*C)/(a^2 + b^2) - c*z + (a^2*c*z)/(a^2 + b^2) - (a*sqrt(-b^2*(C^2 - a^2*S^2 - b^2*S^2 - 2*c*C*z + c^2*z^2)))/(a^2 + b^2))/b

    return Float64(x₁), Float64(y₁), Float64(x₂), Float64(y₂)
end

##

#sphere radius
R = 1.0
#line colatitude
Θ = π/4
#spherical circle representing crater
ϕ = 2π*rand()
θ = Θ + (2*rand() - 1)/3
r = rand()/2

figure()
makegrid()
plotcolat(Θ)
plotcirc(θ, ϕ, r, color=:C0)

if abs(θ - Θ) < r/R
    println("intersection!")
    x₁, y₁, x₂, y₂ = intersection(big(θ), big(ϕ), big(r), big(Θ), big(R))
    z = cos(Θ)
    plotpoint(x₁, y₁, z, color="C3")
    plotpoint(x₂, y₂, z, color="C3")
    println(Θ - cart2usph(x₁, y₁, z)[1])
    println(Θ - cart2usph(x₂, y₂, z)[1])
else
    println("no intersection")
end