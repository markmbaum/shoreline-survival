using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using LinearAlgebra

include(scriptsdir("plotting_functions.jl"))

##

#sphere radius
R = 10.0
#spherical circle representing crater
θ, ϕ = sphrand()
r = rand()*R
#segment which might intersect
s = SphericalSegment(sphrand(), sphrand())

#prepare plot
figure()
makegrid(R)
plotcirc(θ, ϕ, r, R=R, color=:C0)
plotgreatcirc(s, R=R)
plotseg(s, R=R, color=:C1, linewidth=2)
plotpoint(0, 0, 0)

#find intersection (or not)
P0 = sph2cart(θ, ϕ)
P1, P2 = sph2cart(s.a), sph2cart(s.b)
N = unit(P1 × P2)
𝓁 = r/R
ψ = asin(P0 ⋅ N)
if abs(ψ) <= 𝓁
    println("intersection!")
    P2′ = unit(P2 - P1*(P1 ⋅ P2))
    A = P0 ⋅ P1
    B = P0 ⋅ P2′
    C = cos(𝓁)
    t₀ = atan(B,A)
    Δt = acos(cos(𝓁)/sqrt(A^2 + B^2))
    t₁ = t₀ + Δt
    t₂ = t₀ - Δt
    plotvec(cos(t₁)*P1 + sin(t₁)*P2′, R=R, color=:C3)
    plotvec(cos(t₂)*P1 + sin(t₂)*P2′, R=R, color=:C3)
else
    println("no intersection")
end

tight_layout()
