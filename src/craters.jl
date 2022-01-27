#==============================================================================
The following definitions handle crater production 
==============================================================================#

export ♂ᵣ, ♂ₐ
export agescaling, craterdensities, cratercounts

#mean radius of Mars [m]
const ♂ᵣ = 3.3895e6
#surface area of Mars [m^2]
const ♂ₐ = 4π*(♂ᵣ^2)

#crater size-frequency bins, see Table 1 in:
#Michael, G. G. Planetary surface dating from crater size–frequency distribution measurements: Multiple resurfacing episodes and differential isochron fitting. Icarus 226, 885–890 (2013)
# i - index of bin
# D - minimum diameter of bin [km]
# N - crater density [1/km^2/Gyr]
const 𝐂 = Dict(
    "i" => collect(-16:19),
    "D" => Float64[0.00391, 0.00553, 0.00782, 0.0111, 0.01565, 0.0221, 0.0313, 0.0442, 0.06251, 0.0887, 0.125, 0.177, 0.25, 0.354, 0.5, 0.7075, 1, 1.415, 2, 2.83, 4, 5.66, 8.05, 11.32, 16.05, 22.63, 32.05, 45.3, 64.05, 90.6, 128.05, 181.1, 256.05, 362.1, 512.05, 724.1],
    "N" => Float64[4.04e3, 2.33e3, 1.14e3, 4.58e2, 1.91e2, 6.66e1, 2.40e1, 9.44e0, 3.30e0, 1.22e0, 4.37e-1, 1.47e-1, 4.70e-2, 1.38e-2, 4.02e-3, 1.15e-3, 3.08e-4, 1.28e-4, 6.85e-5, 3.67e-5, 1.98e-5, 1.06e-5, 5.68e-6, 3.04e-6, 1.62e-6, 8.71e-7, 4.67e-7, 2.40e-7, 1.12e-7, 5.21e-8, 2.43e-8, 1.13e-8, 5.28e-9, 2.47e-9, 1.15e-9, 5.37e-10]
)

#see equation 3 in:
#Michael, G. G. Planetary surface dating from crater size–frequency distribution measurements: Multiple resurfacing episodes and differential isochron fitting. Icarus 226, 885–890 (2013)
𝒻scale(gya) = 3.79e-14*(exp(6.93*gya) - 1) + 5.84e-4*gya
agescaling(gya) = 𝒻scale(gya)/𝒻scale(1)

meanradius(bin::Int) = 1e3*exp2(bin/2 + 1/4)/2

function craterdensities(gya)
    #mean crater radius for each bin [meters]
    r = meanradius.(𝐂["i"])
    #frequency/density [craters/m^2]
    ρ = agescaling(gya)*𝐂["N"]/1e6
    return r, ρ
end

function cratercounts(gya, area)
    #radius bins and frequencies
    r, ρ = craterdensities(gya)
    #counts [craters]
    n = ρ*area
    return r, ρ, n
end

#--------------------------------------
abstract type Crater{T<:AbstractFloat} end

#multiplcation scales crater radius by a factor of f, returning a new GlobalCrater
function *(c::Crater{T}, f::Real) where {T}
    @assert f >= 0 "crater radius cannot be negative"
    GlobalCrater(c.θ, c.ϕ, T(f*c.r))
end

#--------------------------------------
export GlobalCrater

struct GlobalCrater{T} <: Crater{T}
    θ::T
    ϕ::T
    r::T
end

function Base.show(io::IO, c::GlobalCrater{T}) where {T}
    θ, ϕ, r = map(x->round(x, sigdigits=4), (c.θ, c.ϕ, c.r))
    print(io, "GlobalCrater{$T} θ=$θ, ϕ=$ϕ, r=$r")
end

#creates a randomly located crater with radius r
function randglobalcrater(r::T) where {T<:AbstractFloat}
    θ, ϕ = sphrand()
    GlobalCrater(T(θ), T(ϕ), r)
end

#creates a randomly located crater with radius r, using a specific random number generator
function randglobalcrater(r::T, rng::AbstractRNG) where {T<:AbstractFloat}
    θ, ϕ = sphrand(rng)
    GlobalCrater(T(θ), T(ϕ), r)
end

#computes spherical distance between crater center and coordinate [θ,ϕ] using sphere radius r
sphdist(c::GlobalCrater, θ, ϕ, R=♂ᵣ) = sphdist(c.θ, c.ϕ, θ, ϕ, R)

#draws a circle on the sphere representing crater boundary
sphcirc(c::GlobalCrater, R=♂ᵣ; N::Int=51) = sphcirc(c.θ, c.ϕ, c.r, R; N=N)

#--------------------------------------
export LocalCrater

struct LocalCrater{T} <: Crater{T}
    x::T
    y::T
    r::T
end

function Base.show(io::IO, c::LocalCrater{T}) where {T}
    x, y, r = map(x->round(x, sigdigits=4), (c.x, c.y, c.r))
    print(io, "LocalCrater{$T} x=$x, y=$y, r=$r")
end

#creates a randomly located crater with radius r
function randlocalcrater(r::T, Δx::Real, Δy::Real) where {T<:AbstractFloat}
    x = Δx*rand()
    y = Δy*rand()
    LocalCrater{T}(T(x), T(y), r)
end

#creates a randomly located crater with radius r, using a specific random number generator
function randlocalcrater(r::T, Δx::Real, Δy::Real, rng::AbstractRNG) where {T<:AbstractFloat}
    x = Δx*rand(rng)
    y = Δy*rand(rng)
    LocalCrater{T}(T(x), T(y), r)
end

#----------------------------------------
#general functions for crater populations

abstract type CraterPopulation end

function populatebins(t::Real, area::Real, rmin::Real, nmax::Real)
    @assert t > 0
    @assert rmin >= 0
    @assert nmax >= 0
    #mean crater radius [m] and count for each bin
    r, _, n = cratercounts(t, area)
    #ROUND DOWN, just to be conservative
    counts = n .|> floor .|> Int64
    #don't use bins with huge counts for now (or zero craters)
    idx = @. (counts <= nmax) & (counts > 0) & (r >= rmin)
    r, counts = r[idx], counts[idx]
    #iterate more efficiently with larger craters in the first bins
    idx = sortperm(r, rev=true)
    r = r[idx]
    counts = counts[idx]
    return r, counts
end

function Base.iterate(P::CraterPopulation, state::NTuple{2,Int64}=(1,1))
    @unpack bins, counts, r, rng = P
    #current bin and index within that bin
    bin, idx = state
    #stop at the end of the last bin
    bin > bins && return nothing
    #radius is always the same within a bin for this new crater
    @inbounds c = generatecrater(P, r[bin])
    #continue in the same bin or start the next one
    @inbounds (idx > counts[bin]) ? (c, (bin+1,1)) : (c, (bin,idx+1))
end

function Base.show(io::IO, P::CraterPopulation)
    r = map(r->round(r, sigdigits=8), P.r)
    counts = [format(c, commas=true) for c in P.counts] 
    pretty_table(io,
        Any[1:P.bins r counts],
        ["", "Mean Radius [m]", "Count"]
    )
    N = format(P.N, commas=true)
    print(io, "Total craters: $N\n")
end

#--------------------------------------
export GlobalPopulation

#==============================================================================
This is a very low memory representation of a global crater population.
Crater coordinates are not stored, but produced on the fly when iterating.
A GlobalPopulation behaves like a random crater generator, so each
will have it's own random number generator for seeding and
reproducibility. Individual craters can only be generated via iteration,
not by indexing or anything else
===============================================================================#
struct GlobalPopulation <: CraterPopulation
    bins::Int64 #number of radius bins
    counts::Vector{Int64} #crater count in each bin
    N::Int64 #total number of craters
    r::Vector{Float64} #mean radius of each bin
    rng::Xoshiro
end

function GlobalPopulation(r::Vector, counts::Vector{Int64}, seed=1)
    @assert length(r) == length(counts)
    GlobalPopulation(length(r), counts, sum(counts), r, Xoshiro(seed))
end

function GlobalPopulation(t::Real; rmin::Real=0, nmax::Real=Inf, seed=1)
    #generate crater size-count bins
    r, counts = populatebins(t, ♂ₐ, rmin, nmax)
    #store it all in, effectively, a random crater generator
    return GlobalPopulation(r, counts, seed)
end

Base.eltype(::GlobalPopulation) = GlobalCrater{Float64}

Base.length(P::GlobalPopulation) = P.N

generatecrater(P::GlobalPopulation, r) = randglobalcrater(r, P.rng)

#--------------------------------------
export LocalPopulation

#==============================================================================
This is a similarly low memory representation of a crater population, but for 
a limited cartesian domain instead of a global spherical domain.
===============================================================================#

struct LocalPopulation <: CraterPopulation
    bins::Int64 #number of radius bins
    counts::Vector{Int64} #crater count in each bin
    N::Int64 #total number of craters
    r::Vector{Float64} #mean radius of each bin
    rng::Xoshiro #random number generator
    Δx::Float64 #horizontal size of domain [meters]
    Δy::Float64 #vertical size of domain [meters]
end

function LocalPopulation(r::Vector, counts::Vector{Int64}, Δx::Real, Δy::Real, seed=1)
    @assert length(r) == length(counts)
    LocalPopulation(length(r), counts, sum(counts), r, Xoshiro(seed), Δx, Δy)
end

function LocalPopulation(t::Real, Δx::Real, Δy::Real; rmin::Real=0, nmax::Real=Inf, seed=1)
    #generate crater size-count bins
    r, counts = populatebins(t, Δx*Δy, rmin, nmax)
    #store it all in, effectively, a random crater generator
    return LocalPopulation(r, counts, Δx, Δy, seed)
end

Base.eltype(::LocalPopulation) = LocalCrater{Float64}

Base.length(P::LocalPopulation) = P.N

generatecrater(P::LocalPopulation, r) = randlocalcrater(r, P.Δx, P.Δy, P.rng)
