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
rₑ = 1.5
Δ = 50
rmin = 100
nmax = Inf

##

N = 21
seglens = Vector{Vector{Float64}}(undef, N)
gaplens = Vector{Vector{Float64}}(undef, N)
@threads for i ∈ 1:N
    res = globalsimulation(t, θₛ, rₑ, Δ, rmin=rmin, nmax=nmax, seed=i)
    seglens[i] = segdistances(res, θₛ)
    gaplens[i] = gapdistances(res, θₛ)
end
seglen = vcat(seglens...)
gaplen = vcat(gaplens...);

##

println("percentiles")
for q ∈ [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
    s = quantile(seglen, q)
    g = quantile(gaplen, q)
    n = Int64(100*q)
    println("  $(n)\n    seg: $s\n    gap: $g")
end

##

fig = figure(figsize=(8,3))

subplot(1,2,1)
hist(seglen/1e3, bins=50, log=true, alpha=0.6, color=:black)
xlabel("Segment Length [km]")
ylabel("Count")

subplot(1,2,2)
hist(gaplen/1e3, bins=50, log=true, alpha=0.6, color=:black)
xlabel("Gap Length [km]")

tight_layout()
fig[:savefig](plotsdir("results", "histograms"), dpi=500)
