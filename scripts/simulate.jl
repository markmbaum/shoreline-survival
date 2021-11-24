using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using IterTools: product
using Base.Threads: @threads
using Statistics

## functions

function batch(t, θₛ, rₑ, rmin, nmax, N)::Vector{SimulationResult}
    res = Vector{SimulationResult}(undef, N)
    @threads for i = 1:N
        res[i] = simulateimpacts(t, θₛ, rₑ, rmin=rmin, nmax=nmax, seed=i)
    end
    return res
end

function simulate(params, N::Int, rmin, nmax, fn::String)::Nothing

    #write column names to file
    open(fn, "w") do io
        println(io, "t,theta,re,f,segmean,segmedian,segmax")
    end
    
    #do simulations in parallel batches, writing to file along the way
    for (t, θₛ, rₑ) ∈ params
        #run many simulations in parallel
        results = batch(t, θₛ, rₑ, rmin, nmax, N)
        #append the results to file
        open(fn, "a") do io
            for result ∈ results
                f = result.destroyed
                seglen = segmentlengths(result)
                segmean = mean(seglen)
                segmedian = median(seglen)
                segmax = maximum(seglen)
                println(io, "$t,$θₛ,$rₑ,$f,$segmean,$segmedian,$segmax")
            end
        end
    end
    return nothing
end

## parameter selection/definition

#times [Gya]
t = [LinRange(4, 3.5, 21); LinRange(3.45, 3, 10)]
#shoreline colatitudes [rad]
θₛ = [π/5, π/4, π/3, π/2]
#ejecta distance as multiple of radius
rₑ = [1.0, 1.25, 1.5, 1.75, 2.0]
#minimum crater radius
rmin = 250
#maximum number of craters per bin (should be a HIGH ceiling)
nmax = Inf

#create parameter combinations
params = collect(product(t, θₛ, rₑ));

##

#number of trials for each parameter combination
N = 25000
#simulate and write to file all at once
simulate(params, N, rmin, nmax, datadir("sims", "simulations.csv"))

