using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using IterTools: product
using Base.Threads: @threads
using Statistics

## functions

function batch(t, segments, rₑ, Δ, rmin, nmax, N)::Vector{SimulationResult{SphericalSegment}}
    res = Vector{SimulationResult{SphericalSegment}}(undef, N)
    @threads for i = 1:N
        res[i] = simulateimpacts(t, segments, rₑ, Δ, rmin=rmin, nmax=nmax, seed=i)
    end
    return res
end

sigdig(x) = round(x, sigdigits=6)

function simulate(params, segments, N::Int, rmin, nmax, fn::String)::Nothing

    #write column names to file
    open(fn, "w") do io
        println(io, "t,re,overlap,f,impacts,segmean,segmedian,segmax")
    end
    
    #do simulations in parallel batches, writing to file along the way
    count = 1
    L = length(params)
    for (t, rₑ, Δ) ∈ params
        #run many simulations in parallel
        results = batch(t, segments, rₑ, Δ, rmin, nmax, N)
        #print a little update
        println(stdout, "batch $count/$L complete")
        flush(stdout)
        count += 1
        #append results to the csv file
        open(fn, "a") do io
            for result ∈ results
                d = segmentdistances(result)
                print(io,
                    sigdig(t), ',',
                    sigdig(rₑ), ',',
                    sigdig(Δ), ',',
                    sigdig(result.destroyed), ',',
                    result.impacts, ',',
                    sigdig(mean(d)), ',',
                    sigdig(median(d)), ',',
                    sigdig(maximum(d)), '\n'
                )
            end
        end
    end
    return nothing
end

## parameter selection/definition

#times [Gya], denser at older periods
t = [LinRange(4, 3.75, 11); LinRange(3.7, 3.5, 5); LinRange(3.4, 3, 5)]
#ejecta distance as multiple of radius [-]
rₑ = [1.0, 1.5, 2.0]
#required overlap distance [m]
Δ = [5.0] #[5e0, 5e1, 5e2]
#minimum crater radius [m]
rmin = 200
#maximum number of craters per bin (should be a HIGH ceiling)
nmax = Inf

#create parameter combinations
params = collect(product(t, rₑ, Δ));
println("$(length(params)) parameter combinations")

## load the putative shoreline

segments = readsegments(
    datadir("exp_pro", "parker_1989_contact_1a.csv"),
    minarc=0.05
);
println("$(length(segments)) initial segments")

##
 
#simulate and write to file all at once
simulate(
    params,
    segments,
    48, #should be at least the number of available threads
    rmin,
    nmax,
    datadir("sims", "mapped.csv")
)
