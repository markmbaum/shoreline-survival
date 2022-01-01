using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using Statistics: quantile
using PyPlot
using Base.Threads: @threads

pygui(true)

##

t = 4
θₛ = π/3
rₑ = 1.0
Δ = 50
rmin=500
nmax=Inf

##

#shoreline coordinates
#fn = datadir("exp_pro", "parker_1989_contact_1a.csv")
#read the coordinates into segments with appropriate spacing
#segments = readsegments(fn, minarc=0.01);

##

seglens = Vector{Float64}[]
@threads for seed ∈ 1:100
    res = simulateimpacts(t, θₛ, rₑ, Δ, rmin=rmin, nmax=nmax, seed=seed)
    push!(seglens, segmentdistances(res, θₛ))
end
seglen = vcat(seglens...);

##

println("50th percentile: ", quantile(seglen, 0.5))
println("75th percentile: ", quantile(seglen, 0.75))
println("95th percentile: ", quantile(seglen, 0.95))
println("99th percentile: ", quantile(seglen, 0.99))

##

figure(figsize=(8,3))
subplot(1,2,1)
hist(seglen/1e3, bins=50, log=true, density=false, alpha=0.6, color=:black)
xlabel("Segment Length [km]")
ylabel("Count");
subplot(1,2,2)
hist(seglen/1e3, bins=50, log=false, density=false, alpha=0.6, color=:black)
xlabel("Segment Length [km]")
ylabel(nothing);
tight_layout()

##

figure(figsize=(4,3))
hist(seglen/1e3, bins=50, log=true, density=true, alpha=0.6, color=:black)
xlabel("Segment Length [km]")
ylabel("Density");
tight_layout()

figure(figsize=(4,3))
hist(seglen/1e3, bins=50, log=false, density=true, alpha=0.6, color=:black)
xlabel("Segment Length [km]")
ylabel("Density");
tight_layout()