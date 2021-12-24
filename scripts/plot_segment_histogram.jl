using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using PyPlot
using Base.Threads: @threads

pygui(true)

##

t = 3
θₛ = π/3
rₑ = 1.5
Δ = 50
rmin=100
nmax=1e8

##

seglens = Vector{Float64}[]
@threads for seed ∈ 1:100
    res = simulateimpacts(t, θₛ, rₑ, Δ, rmin=rmin, nmax=nmax, seed=seed)
    push!(seglens, segmentdistances(res, θₛ))
end
seglen = vcat(seglens...);

##

figure()
subplot(1,2,1)
hist(seglen/1e3, bins=50, log=true, density=false, alpha=0.6, color=:black)
xlabel("Segment Length [km]")
ylabel("Count");
subplot(1,2,2)
hist(seglen/1e3, bins=50, log=false, density=false, alpha=0.6, color=:black)
xlabel("Segment Length [km]")
ylabel(nothing);
tight_layout()