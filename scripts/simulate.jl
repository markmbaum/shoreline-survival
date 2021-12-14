using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using IterTools: product
using Base.Threads: @threads
using Statistics

## functions

function batch(t, θₛ, rₑ, Δ, rmin, nmax, N)::Vector{SimulationResult}
    res = Vector{SimulationResult}(undef, N)
    @threads for i = 1:N
        res[i] = simulateimpacts(t, θₛ, rₑ, Δ, rmin=rmin, nmax=nmax, seed=i)
    end
    return res
end

sigdig(x) = round(x, sigdigits=6)

function simulate(params, N::Int, rmin, nmax, fn::String)::Nothing

    #write column names to file
    open(fn, "w") do io
        println(io, "t,theta,re,overlap,f,segmean,segmedian,segmax")
    end
    
    #do simulations in parallel batches, writing to file along the way
    count = 1
    for (t, θₛ, rₑ, Δ) ∈ params
        #run many simulations in parallel
        results = batch(t, θₛ, rₑ, Δ, rmin, nmax, N)
        #print a little update
        println(stdout, "batch $count/$(length(params)) complete")
        flush(stdout)
        count += 1
        #append results to the csv file
        open(fn, "a") do io
            for result ∈ results
                seglen = segmentlengths(result, θₛ)
                print(io,
                    sigdig(t), ',',
                    sigdig(θₛ), ',',
                    sigdig(rₑ), ',',
                    sigdig(Δ), ',',
                    sigdig(result.destroyed), ',',
                    sigdig(mean(seglen)), ',',
                    sigdig(median(seglen)), ',',
                    sigdig(maximum(seglen)), '\n'
                )
            end
        end
    end
    return nothing
end

## parameter selection/definition

#times [Gya], denser at older periods
t = [LinRange(4, 3.75, 11); LinRange(3.7, 3.5, 5); LinRange(3.4, 3, 5)]
#shoreline colatitudes [rad]
θₛ = [π/5, π/4, π/3, π/2]
#ejecta distance as multiple of radius
rₑ = [1.0, 1.25, 1.5, 1.75, 2.0]
#required overlap distance
Δ = [5e0, 5e1, 5e2]
#minimum crater radius
rmin = 250
#maximum number of craters per bin (should be a HIGH ceiling)
nmax = Inf

#create parameter combinations
params = collect(product(t, θₛ, rₑ, Δ));

##
 
#simulate and write to file all at once
simulate(
    params,
    10000, #number of trials for each parameter combination
    rmin,
    nmax,
    datadir("sims", "simulations.csv")
)
