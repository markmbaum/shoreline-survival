using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using PyPlot

pygui(true)

##

function plotvec(θ, ϕ, color="C0")
    x, y, z = sph2cart(θ, ϕ)
    plot3D([0,x], [0,y], [0,z], color=color)
end

plotvec(p::SphericalPoint, color="C0") = plotvec(p.θ, p.ϕ, color)

function plotpoint(θ, ϕ, color="C0")
    x, y, z = sph2cart(θ, ϕ)
    plot3D([x], [y], [z], ".", color=color)
end

plotpoint(p::SphericalPoint, color="C0") = plotpoint(p.θ, p.ϕ, color)

function segplot(s::SphericalSegment, color="C0")
    θ = LinRange(s.a.θ, s.b.θ, 10)
    ϕ = LinRange(s.a.ϕ, s.b.ϕ, 10)
    x, y, z = sph2cart(θ, ϕ, 1.0)
    plot3D(x, y, z, color=color)
end

function craterplot(c::Crater, color="C0")
    x, y, z = sphcirc(c, 1)
    plot3D(x, y, z, color=color)
end

function setlim()
    xlim(-1.05,1.05)
    ylim(-1.05,1.05)
    zlim(-1.05,1.05)
end

function makegrid(N::Int=6, L::Int=100)
    ϕ = LinRange(0, 2π, L)
    for θ in LinRange(0, π, N+2)[2:end-1]
        r = sin(θ)
        x = r*sin.(ϕ)
        y = r*cos.(ϕ)
        z = fill(cos(θ), L)
        plot3D(x, y, z, "k", linewidth=0.7, alpha=0.2)
    end
    θ = LinRange(0, π, L ÷ 2)
    for ϕ in LinRange(0, 2π, 2N)[1:end-1]
        x, y, z = sph2cart(θ, fill(ϕ, L ÷ 2), 1.0)
        plot3D(x, y, z, "k", linewidth=0.7, alpha=0.2)
    end
    setlim()
end

##

figure()

a = SphericalPoint(π/4, 1)
b = SphericalPoint(3π/4, 1.5)
plotpoint(a)
plotpoint(b, "C3")
x = sphrand()
plotpoint(x[1], x[2], "C1")
y = sphrand()
plotpoint(y[1], y[2], "C2")

N = 50
θ = LinRange(a.θ, b.θ, N)
ϕ = LinRange(a.ϕ, b.ϕ, N)
for i = 1:N
    plotvec(θ[i], ϕ[i])
    z = pointrotation(x[1], x[2], a, SphericalPoint(θ[i], ϕ[i]))
    plotvec(z[1], z[2], "C1")
    z = pointrotation(y[1], y[2], a, SphericalPoint(θ[i], ϕ[i]))
    plotvec(z[1], z[2], "C2")
end

z = pointrotation(x[1], x[2], a, b)
plotpoint(z[1], z[2], "C3")
z = pointrotation(y[1], y[2], a, b)
plotpoint(z[1], z[2], "C3")

makegrid()

##

figure()

c = Crater(0.05)
θ, ϕ = sphrand()
p = SphericalPoint(c.θ, c.ϕ)
s = SphericalSegment(θ, ϕ, θ + rand()/4, ϕ + rand()/4)
craterplot(c)
segplot(s)
println("arclength 1: $(arclength(s.a, p))")
println("arclength 2: $(arclength(s.b, p))")

craterplot(Crater(0.0, 0.0, c.r), "C1")
r = polerotation(s, c.θ, c.ϕ)
segplot(r, "C1")
println("colatitude 1: $(r.a.θ)")
println("colatitude 2: $(r.b.θ)")

makegrid()
