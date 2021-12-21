using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using BenchmarkTools
using PyPlot
using MultiAssign

pygui(true)

##

function plotcirc(θ, ϕ, r; R=1, color="C0", linewidth=1.5, N=1000)
    x, y, z = sphcirc(θ, ϕ, r, R, N=N)
    plot3D(x, y, z, color=color, linewidth=linewidth)
    nothing
end

function plotpoint(x, y, z; R=1, color="C3", size=5)
    plot3D(x, y, z, ".", color=color, markersize=size)
    nothing
end

function plotpoint(p::SphericalPoint; kw...)
    plotpoint(sph2cart(p); kw...)
end

function plotvec(v, w)
    plot3D(
        [v[1], w[1]],
        [v[2], w[2]],
        [v[3], w[3]]
    )
end

function plotvec(v)
    plot3D(
        [0.0, v[1]],
        [0.0, v[2]],
        [0.0, v[3]]
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

function plotcolat(θ; color="k", linewidth=1, N=250, R=1)
    @multiassign x, y, z = zeros(N)
    ϕ = LinRange(0, 2π, N)
    for i ∈ 1:N
        x[i], y[i], z[i] = sph2cart(θ, ϕ[i])
    end
    plot3D(R*x, R*y, R*z, color=color, linewidth=linewidth)
    nothing
end

function setlim()
    xlim(-1.05,1.05)
    ylim(-1.05,1.05)
    zlim(-1.05,1.05)
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
    setlim()
    axis("off")
    nothing
end

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
r = rand()

figure()
makegrid()
plotcolat(Θ)
plotcirc(θ, ϕ, r)

if abs(θ - Θ) < r/R
    println("intersection!")
    x₁, y₁, x₂, y₂ = intersection(big(θ), big(ϕ), big(r), big(Θ), big(R))
    z = cos(Θ)
    plotpoint(x₁, y₁, z)
    plotpoint(x₂, y₂, z)
    println(Θ - cart2usph(x₁, y₁, z)[1])
    println(Θ - cart2usph(x₂, y₂, z)[1])
else
    println("no intersection")
end