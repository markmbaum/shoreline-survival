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

##

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

xlim(-1.1,1.1)
ylim(-1.1,1.1)
zlim(-1.1,1.1)