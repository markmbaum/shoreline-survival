using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using LinearAlgebra
using MultiAssign
using PyPlot

pygui(true)

##

function plotcirc(θ, ϕ, r; R=1, color="C0", linewidth=1.5, N=1000)
    x, y, z = sphcirc(θ, ϕ, r, R, N=N)
    plot3D(x, y, z, color=color, linewidth=linewidth)
    nothing
end

function plotpoint(p::SphericalPoint; R=1, color="C3", size=5)
    x, y, z = sph2cart(p)
    plot3D(R*x, R*y, R*z, ".", color=color, markersize=size)
    nothing
end

function plotvec(v, w; R=1)
    plot3D(
        R*[v[1], w[1]],
        R*[v[2], w[2]],
        R*[v[3], w[3]]
    )
end

function plotvec(v; R=1)
    plot3D(
        R*[0.0, v[1]],
        R*[0.0, v[2]],
        R*[0.0, v[3]]
    )
end

function plotseg(s::SphericalSegment; R=1, N=1000, color="C3", linewidth=1.5)
    plotpoint(s.a, R=R, color=color)
    plotpoint(s.b, R=R, color=color)
    C = GreatCircle(s)
    t = LinRange(0.0, arclength(s), N)
    @multiassign x, y, z = zeros(N)
    for i ∈ 1:N
        x[i], y[i], z[i] = C(t[i])
    end
    plot3D(R*x, R*y, R*z, color=color, linewidth=linewidth)
    nothing
end

function plotgreatcirc(s::SphericalSegment; R=1, N=1000, color="k", linewidth=0.75)
    C = GreatCircle(s)
    t = LinRange(0.0, 2π, N)
    @multiassign x, y, z = zeros(N)
    for i ∈ 1:N
        x[i], y[i], z[i] = C(t[i])
    end
    plot3D(R*x, R*y, R*z, color=color, linewidth=linewidth)
    nothing
end

function setlim(R)
    xlim(-1.05*R, 1.05*R)
    ylim(-1.05*R, 1.05*R)
    zlim(-1.05*R, 1.05*R)
    nothing
end

function makegrid(R=1.0, N::Int=6, L::Int=100)
    ϕ = LinRange(0, 2π, L)
    for θ in LinRange(0, π, N+2)[2:end-1]
        r = sin(θ)
        x = r*sin.(ϕ)
        y = r*cos.(ϕ)
        z = fill(cos(θ), L)
        plot3D(R*x, R*y, R*z, "k", linewidth=0.7, alpha=0.2)
    end
    θ = LinRange(0, π, L ÷ 2)
    for ϕ in LinRange(0, 2π, 2N)[1:end-1]
        x, y, z = sph2cart(θ, fill(ϕ, L ÷ 2), 1.0)
        plot3D(R*x, R*y, R*z, "k", linewidth=0.7, alpha=0.2)
    end
    setlim(R)
    axis("off")
    nothing
end

##

#sphere radius
R = 10.0
#spherical circle representing crater
θ, ϕ = sphrand()
r = rand()*R
#segment which might intersect
s = SphericalSegment(sphrand(), sphrand())



figure()
makegrid(R)
plotcirc(θ, ϕ, r, R=R)
plotgreatcirc(s, R=R)
plotseg(s, R=R)

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
    plotvec(cos(t₁)*P1 + sin(t₁)*P2′, R=R)
    plotvec(cos(t₂)*P1 + sin(t₂)*P2′, R=R)
else
    println("no intersection")
end

tight_layout()
