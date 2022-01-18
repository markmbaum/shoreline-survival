using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using IterTools: product
using Statistics

##-----------------------------------------------------------------------------
# functions

sigdig(x) = round(x, sigdigits=6)

function simulate(params, segments, N::Int, rmin, nmax, fn::String)::Nothing

    #write column names to file
    open(fn, "w") do io
        println(io,
            "seed,",
            "t,",
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
    count = 1
    L = length(params)
    for (t, rₑ, Δ) ∈ params
        #trials/realizations
        println(stdout, "batch $count/$L beginning")
        flush(stdout)
        count += 1
        for i ∈ 1:N
            #run a simulation
            result = globalsimulation(t, segments, rₑ, Δ, rmin=rmin, nmax=nmax, seed=i)
            println("  trial $i complete")
            flush(stdout)
            #append results to the csv file
            open(fn, "a") do io
                ds = segdistances(result)
                dg = gapdistances(result)
                print(io,
                    i, ',',
                    sigdig(t), ',',
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
t = [LinRange(4, 3.6, 5); 3.3; 3.0]
#ejecta distance as multiple of radius [-]
rₑ = [1.0, 1.5, 2.0]
#required overlap distance [m]
Δ = [50.0]
#minimum crater radius [m]
rmin = 100
#maximum number of craters per bin (should be a HIGH ceiling)
nmax = Inf
#number of simulations for each parameter combo
N = 10 #doesn't have to be a multiple of thread count here

##-----------------------------------------------------------------------------
# MAIN

#create parameter combinations
params = collect(product(t, rₑ, Δ));
println(stdout, "$(length(params)) parameter combinations")

## load the putative shoreline

segments = readsegments(
    datadir("exp_pro", "parker_1989_contact_1a.csv"),
    minarc=0.04
);
println(stdout, "$(length(segments)) initial segments")

##

#simulate and write to file all at once
simulate(
    params,
    segments,
    N,
    rmin,
    nmax,
    datadir("sims", "mapped.csv")
)
