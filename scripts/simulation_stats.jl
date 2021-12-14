using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using DataFrames
using CSV
using IterTools: product
using Statistics: quantile
using Base.Threads: @threads

## load the big table of results

sims = CSV.read(datadir("sims", "simulations.csv"), DataFrame)

## get unique values for each parameter

uparams(df, param) = sort(unique(df[!,param]))

t = uparams(sims, :t)
θ = uparams(sims, :theta)
r = uparams(sims, :re)
Δ = uparams(sims, :overlap)

## start a new dataframe for output statistics

qcol(col, q) = col*string(Int(100*q))

params = collect(product(t, θ, r, Δ))
N = length(params)
stats = DataFrame()
for col ∈ [:t, :θ, :r, :Δ]
    stats[!,col] = zeros(N)
end

#columns to take quantiles for
cols = ["f", "segmean", "segmedian", "segmax"]
#quantiles to compute
quan = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]

for col ∈ cols, q ∈ quan
    stats[!,qcol(col,q)] = zeros(N)
end

## compute stats for each parameter combination

@threads for i ∈ 1:N
    #unpack parameter values
    t, θ, r, Δ = params[i]
    #fill parameter values
    stats[i,:t] = t
    stats[i,:θ] = θ
    stats[i,:r] = r
    stats[i,:Δ] = Δ
    #filter simulations
    sl = sims[(sims.t .≈ t) .& (sims.theta .≈ θ) .& (sims.re .≈ r) .& (sims.overlap .≈ Δ), :]
    #compute quantiles for each results column
    for col ∈ cols
        x = quantile(sl[!,col], quan)
        for (j,q) ∈ enumerate(quan)
            stats[i,qcol(col,q)] = x[j]
        end
    end
end

## write results to file

CSV.write(datadir("exp_pro", "simulation_stats.csv"), stats)
