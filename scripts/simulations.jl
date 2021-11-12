using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using DataFrames
using IterTools: product
using Base.Threads: @spawn, Task, fetch
using Statistics
using ProgressMeter
using CSV

##

q2col(q, c) = c*string(Int(100*q))

function simulate(t, θₛ, rmin, nmax, rₑ, N)
    tasks = Vector{Task}(undef,N)
    for i ∈ 1:N
        tasks[i] = @spawn simulateimpacts(t, θₛ, rmin=rmin, nmax=nmax, rₑ=rₑ, seed=i)
    end
    SimulationResult[fetch(task) for task ∈ tasks]
end

## parameter selection/definition

#times [Gya]
t = [LinRange(4, 3.5, 21); LinRange(3.45, 3, 10)]
#shoreline colatitudes [rad]
θₛ = [π/5, π/4, π/3, π/2]
#minimum crater radius
rmin = 250
#maximum number of craters per bin (should be a HIGH ceiling)
nmax = Inf
#ejecta distance as multiple of radius
rₑ = [1.0, 1.25, 1.5, 1.75, 2.0]

#create parameter combinations
params = collect(product(t, θₛ, rₑ));

##

#number of draws per parameter combo
N = 250
#length of parameter combos
L = length(params)
#quantiles to compute
Q = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]

##

#initialize a dataframe with columns for fraction destroyed
df = DataFrame()
for col ∈ [:t, :theta, :re, :fmean, :fstd, :smaxmean, :smaxstd, :smeanmean, :smeanstd]
    df[!,col] = zeros(L)
end
#add columns for quantiles of destroyed fraction and segment length max & mean
for i ∈ 1:length(Q)
    df[!,q2col(Q[i], 'f')] = zeros(L)
    df[!,q2col(Q[i], "smax")] = zeros(L)
    df[!,q2col(Q[i], "smean")] = zeros(L)
end

##

#simulate!
@showprogress 3 for (i, (t, θₛ, rₑ)) ∈ enumerate(params)

    #put parameters in the df
    df[i,:t] = t
    df[i,:theta] = θₛ
    df[i,:re] = rₑ

    #run many simulations in parallel
    results = simulate(t, θₛ, rmin, nmax, rₑ, N)

    #pull out fraction destroyed for each run
    destroyed = map(res->res.destroyed, results)
    #stats on destroyed fraction
    df[i,:fmean] = mean(destroyed)
    df[i,:fstd] = std(destroyed)
    q = quantile(destroyed, Q)
    for j in 1:length(q)
        col = q2col(Q[j], 'f')
        df[i,col] = q[j]
    end

    #get segment lengths of each run
    seglen = map(segmentlengths, results)
    #pull out max and mean segment lengths for each run
    smax = map(maximum, seglen)
    smean = map(mean, seglen)
    #store mean and std of both
    df[i,:smaxmean] = mean(smax)
    df[i,:smaxstd] = std(smax)
    df[i,:smeanmean] = mean(smean)
    df[i,:smeanstd] = std(smean)
    #quantiles
    q = quantile(smax, Q)
    for j ∈ 1:length(q)
        col = q2col(Q[j], "smax")
        df[i,col] = q[j]
    end
    q = quantile(smean, Q)
    for j ∈ 1:length(q)
        col = q2col(Q[j], "smean")
        df[i,col] = q[j]
    end
    
    #write the results occasionally just in case progess is slow
    (i % 10 == 0) && CSV.write(datadir("sims", "simulations.csv"), df)
end

CSV.write(datadir("sims", "simulations.csv"), df)

