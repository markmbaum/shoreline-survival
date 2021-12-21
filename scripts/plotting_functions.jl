using PyPlot
using MultiAssign
pygui(true)

function plotcirc(θ, ϕ, r; R=1, color="k", linewidth=1, N=1000)
    x, y, z = sphcirc(θ, ϕ, r, R, N=N)
    plot3D(x, y, z, color=color, linewidth=linewidth)
    nothing
end

function plotpoint(x, y, z; R=1, color="k", markersize=5)
    plot3D(R*x, R*y, R*z, ".", color=color, markersize=markersize)
    nothing
end

plotpoint(p::AbstractVector; kw...) = plotpoint(p...; kw...)

plotpoint(p::SphericalPoint; kw...) = plotpoint(sph2cart(p); kw...)

function plotvec(v, w; R=1, color="k", linewidth=1)
    plot3D(
        R*[v[1], w[1]],
        R*[v[2], w[2]],
        R*[v[3], w[3]],
        color=color,
        linewidth=linewidth
    )
end

plotvec(v; kw...) = plotvec(zeros(3), v; kw...)

function plotseg(c::CartesianSegment; R=1, N=100, color="k", linewidth=1)
    C = GreatCircle(c)
    t = LinRange(0.0, arclength(s), N)
    @multiassign x, y, z = zeros(N)
    for i ∈ 1:N
        x[i], y[i], z[i] = C(t[i])
    end
    plot3D(R*x, R*y, R*z, color=color, linewidth=linewidth)
    nothing
end

plotseg(s::SphericalSegment; kw...) = plotseg(CartesianSegment(s); kw...)

function plotgreatcirc(C::GreatCircle; R=1, N=250, color="gray", linewidth=1)
    t = LinRange(0.0, 2π, N)
    @multiassign x, y, z = zeros(N)
    for i ∈ 1:N
        x[i], y[i], z[i] = C(t[i])
    end
    plot3D(R*x, R*y, R*z, color=color, linewidth=linewidth, zorder=-1)
    nothing
end

plotgreatcirc(s::SphericalSegment; kw...) = plotgreatcirc(GreatCircle(s); kw...)

function setlim(R)
    xlim(-1.05*R, 1.05*R)
    ylim(-1.05*R, 1.05*R)
    zlim(-1.05*R, 1.05*R)
    nothing
end

function makegrid(R=1.0, N::Int=5, L::Int=100)
    ϕ = LinRange(0, 2π, L)
    for θ in LinRange(0, π, N+2)[2:end-1]
        r = sin(θ)
        x = r*sin.(ϕ)
        y = r*cos.(ϕ)
        z = fill(cos(θ), L)
        plot3D(R*x, R*y, R*z, "k", linewidth=0.5, alpha=0.15, zorder=-100)
    end
    θ = LinRange(0, π, L ÷ 2)
    for ϕ in LinRange(0, 2π, 2N)[1:end-1]
        x, y, z = sph2cart(θ, fill(ϕ, L ÷ 2), 1.0)
        plot3D(R*x, R*y, R*z, "k", linewidth=0.5, alpha=0.15, zorder=-100)
    end
    setlim(R)
    axis("off")
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

println("plotting functions loaded")