using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using PyPlot

pygui(true)

##

function plotsimulation(t, θₛ, rmin, nmax, seed)
    figure()
    R = ♂ᵣ
    #plot the shoreline ring
    θ = fill(θₛ, 100)
    ϕ = LinRange(0, π, 100)
    x, y, z = sph2cart(θ, ϕ, R)
    plot3D(x, y, z, color="C0", linewidth=0.75)
    #plot all the craters
    P = GlobalPopulation(t, rmin=rmin, nmax=nmax, seed=seed)
    θ = fill(θₛ, 10)
    for crater ∈ P
        #take only one hemisphere
        if 0.05 < crater.ϕ < π - 0.05
            x, y, z = sphcirc(crater)
            if ♂ᵣ*abs(θₛ - crater.θ) < crater.r
                plot3D(x, y, z, color="k", linewidth=0.75)
                ϕ₁, ϕ₂ = intersection(crater, θₛ, ♂ᵣ)
                ϕ₂ += (ϕ₁ > ϕ₂) ? 2π : 0
                ϕ = LinRange(ϕ₁, ϕ₂, 10)
                x, y, z = sph2cart(θ, ϕ, ♂ᵣ)
                plot3D(x, y, z, color="C3", linewidth=1.25, solid_capstyle="round")
            else
                plot3D(x, y, z, color="gray", linewidth=0.5, alpha=0.8)
            end
        end
    end
    axis("off")
    gca()[:view_init](0,90)
end

#these parameters produce some very large craters and few small ones
plotsimulation(4.2, 2π/5, 1000, 1000, 5)

##

savefig(plotsdir("example_simulation.png"), dpi=1200)