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
gaplens = Vector{Float64}[]
@threads for seed ∈ 1:100
    res = simulateimpacts(t, θₛ, rₑ, Δ, rmin=rmin, nmax=nmax, seed=seed)
    push!(seglens, segdistances(res, θₛ))
    push!(gaplens, gapdistances(res, θₛ))
end
seglen = vcat(seglens...)
gaplen = vcat(gaplens...);

##

println("percentiles")
for q ∈ [0.5, 0.75, 0.95, 0.99]
    s = quantile(seglen, q)
    g = quantile(gaplen, q)
    n = Int64(100*q)
    println("  $(n)th\n    seg: $s\n    gap: $g")
end

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