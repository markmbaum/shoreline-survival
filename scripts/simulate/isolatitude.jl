using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using IterTools: product
using Base.Threads: @threads
using Statistics

##-----------------------------------------------------------------------------
# functions

function batch(t, θₛ, rₑ, Δ, rmin, nmax, N)::Vector{SimulationResult{NTuple{2,Float64}}}
    res = Vector{SimulationResult{NTuple{2,Float64}}}(undef, N)
    @threads for i = 1:N
        res[i] = simulateimpacts(t, θₛ, rₑ, Δ, rmin=rmin, nmax=nmax, seed=i)
    end
    return res
end

sigdig(x) = round(x, sigdigits=6)

function simulate(params, N::Int, rmin, nmax, fn::String)::Nothing

    #write column names to file
    open(fn, "w") do io
        println(io,
            "seed,",
            "t,",
            "theta,",
            "re,",
            "overlap,",
            "survived,",
            "destroyed,",
            "impacts,",
            "segmean,",
            "segmedian,",
            "segmax,",
            "segmin,",
            "gapmean,",
            "gapmedian,",
            "gapmax,",
            "gapmin"
        )
    end
    
    #do simulations in parallel batches, writing to file along the way
    flush(stdout)
    count = 1
    L = length(params)
    for (t, θₛ, rₑ, Δ) ∈ params
        #run many simulations in parallel
        results = batch(t, θₛ, rₑ, Δ, rmin, nmax, N)
        #print a little update
        println(stdout, "batch $count/$L complete")
        flush(stdout)
        count += 1
        #append results to the csv file
        open(fn, "a") do io
            for (seed,result) ∈ enumerate(results)
                ds = segdistances(result, θₛ)
                dg = gapdistances(result, θₛ)
                print(io,
                    seed, ',',
                    sigdig(t), ',',
                    sigdig(θₛ), ',',
                    sigdig(rₑ), ',',
                    sigdig(Δ), ',',
                    sigdig(survived(result)), ',',
                    sigdig(destroyed(result)), ',',
                    result.impacts, ',',
                    sigdig(mean(ds)), ',',
                    sigdig(median(ds)), ',',
                    sigdig(maximum(ds)), ',',
                    sigdig(minimum(ds)), ',',
                    sigdig(mean(dg)), ',',
                    sigdig(median(dg)), ',',
                    sigdig(maximum(dg)), ',',
                    sigdig(minimum(dg)), '\n',
                )
            end
        end
    end
    return nothing
end

##-----------------------------------------------------------------------------
# parameter selection/definition

#times [Gya], denser at older periods
t = [LinRange(4, 3.75, 11); LinRange(3.725, 3.5, 5); LinRange(3.4, 3, 5)]
#shoreline colatitudes [rad]
θₛ = map(i->π/i, 2:6)
#ejecta distance as multiple of radius [-]
rₑ = LinRange(1, 2, 9)
#required overlap distance [m]
Δ = [50.0]
#minimum crater radius [m]
rmin = 100
#maximum number of craters per bin (should be a HIGH ceiling)
nmax = Inf
#number of simulations for each parameter combo
N = 144 #should be a multiple of number of available threads

##-----------------------------------------------------------------------------
# MAIN

#create parameter combinations
params = collect(product(t, θₛ, rₑ, Δ));
println(stdout, "$(length(params)) parameter combinations")

##
 
#simulate and write to file all at once
simulate(
    params,
    N,
    rmin,
    nmax,
    datadir("sims", "isolatitude.csv")
)
