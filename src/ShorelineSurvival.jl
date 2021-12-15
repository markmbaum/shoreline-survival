module ShorelineSurvival

using LinearAlgebra: ⋅, ×
using PrettyTables
using UnPack
using Random: AbstractRNG, Xoshiro, seed!
using StaticArrays
using MultiAssign
using Formatting
using Statistics
using Graphs
using SparseArrays
import Base.*

const 𝛕 = 2π

#==============================================================================
This section contains functions for doing various things in spherical geometry.
I use a colatitude coordinate θ (theta) ∈ [0,π] and a longitude
coordinate ϕ (phi) ∈ [0,2π].
==============================================================================#

export sphrand, sph2cart, cart2sph, sphdist, sphcirc

#--------------------------------------

function sphrand(rng::AbstractRNG)::NTuple{2,Float64}
    θ = acos(1 - 2*rand(rng))
    ϕ = 𝛕*rand(rng)
    return θ, ϕ
end

function sphrand()::NTuple{2,Float64}
    θ = acos(1.0 - 2.0*rand())
    ϕ = 𝛕*rand()
    return θ, ϕ
end

function sphrand(n::Int)::NTuple{2,Vector{Float64}}
    @multiassign θ, ϕ = zeros(Float64, n)
    @inbounds for i ∈ 1:n
        θ[i], ϕ[i] = sphrand()
    end
    return θ, ϕ
end

#--------------------------------------

#assumes radius is 1
function sph2cart(θ::T, ϕ::T) where {T<:Real}
    s = sin(θ)
    x = s*cos(ϕ)
    y = s*sin(ϕ)
    z = cos(θ)
    return SVector{3,T}(x, y, z)
end

function sph2cart(θ::T, ϕ::T, r::T) where {T<:Real}
    s = sin(θ)
    x = r*s*cos(ϕ)
    y = r*s*sin(ϕ)
    z = r*cos(θ)
    return SVector{3,T}(x, y, z)
end

function sph2cart(θ::AbstractVector{T}, ϕ::AbstractVector{T}, r::T) where {T}
    @assert length(θ) == length(ϕ)
    @multiassign x, y, z = similar(θ)
    @inbounds for i ∈ 1:length(x)
        x[i], y[i], z[i] = sph2cart(θ[i], ϕ[i], r)
    end
    return x, y, z
end

#--------------------------------------

function cart2sph(x::T, y::T, z::T) where {T<:Real}
    r = sqrt(x*x + y*y + z*z)
    θ = acos(z/r)
    ϕ = (atan(y,x) + 𝛕) % 𝛕
    return SVector{3,T}(θ, ϕ, r)
end

function cart2sph(x::AbstractVector{T},
                  y::AbstractVector{T},
                  z::AbstractVector{T}) where {T}
    @assert length(x) == length(y) == length(z)
    @multiassign θ, ϕ, r = similar(x)
    @inbounds for i ∈ 1:length(θ)
        θ[i], ϕ[i], r[i] = cart2sph(x[i], y[i], z[i])
    end
    return θ, ϕ, r
end

#drops the radius
function cart2usph(x, y, z)
    θ, ϕ, _ = cart2sph(x, y, z)
    return θ, ϕ
end

#--------------------------------------

function arclength(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where {T<:Real}
    v₁ = sph2cart(θ₁, ϕ₁)
    v₂ = sph2cart(θ₂, ϕ₂)
    (v₁ ≈ v₂) ? zero(T) : acos(v₁ ⋅ v₂)
end

sphdist(θ₁, ϕ₁, θ₂, ϕ₂, R) = R*arclength(θ₁, ϕ₁, θ₂, ϕ₂)

function unit(v::SVector{3,T}) where {T}
    x, y, z = v
    L = sqrt(x*x + y*y + z*z)
    return SVector{3,T}(x/L, y/L, z/L)
end

function sphcirc(θ::T, ϕ::T, r::T, R=♂ᵣ; N::Int=75) where {T<:Real}
    #vector from center of sphere to center of circle
    C = sph2cart(θ, ϕ, R)
    #unit vector from sphere center to circle center, normal to circle's plane
    n = unit(C)
    #unit vector perpendicular to n in the x,y plane
    u = SVector{3,T}(-sin(ϕ), cos(ϕ), 0.0)
    #unit vector perpendicular to both n and u using cross product
    v = n × u
    #create vectors of coordinates representing the circle
    @multiassign x, y, z = zeros(T, N)
    u₁, u₂, u₃ = u
    v₁, v₂, v₃ = v
    C₁, C₂, C₃ = C
    @inbounds for (i,ψ) ∈ enumerate(LinRange(0, 𝛕, N))
        s = sin(ψ)
        c = cos(ψ)
        x[i] = C₁ + r*(s*u₁ + c*v₁)
        y[i] = C₂ + r*(s*u₂ + c*v₂)
        z[i] = C₃ + r*(s*u₃ + c*v₃)
    end
    return x, y, z
end

function wrapangle(θ)
    while θ < 0; θ += 𝛕; end
    while θ > 𝛕; θ -= 𝛕; end
    return θ
end

#====
====#

export SphericalPoint, SphericalSegment
export pointrotation

struct SphericalPoint
    θ::Float64
    ϕ::Float64
end

function SphericalPoint(x::NTuple{2,T}) where {T<:Real}
    SphericalPoint(x[1], x[2])
end

function Base.isapprox(a::SphericalPoint, b::SphericalPoint)::Bool
    (a.θ ≈ b.θ) & (a.ϕ ≈ b.ϕ)
end

struct SphericalSegment
    a::SphericalPoint
    b::SphericalPoint
end

sph2cart(p::SphericalPoint)::SVector{3,Float64} = sph2cart(p.θ, p.ϕ)

function arclength(a::SphericalPoint, b::SphericalPoint)::Float64
    arclength(a.θ, a.ϕ, b.θ, b.ϕ)
end

function arclength(s::SphericalSegment)::Float64
    arclength(s.a, s.b)
end

function sphdist(s::SphericalSegment, R::Real=♂ᵣ)::Float64
    R*arclength(s)
end

function pointrotation(θ::Float64, ϕ::Float64, a::SphericalPoint, b::SphericalPoint)::NTuple{2,Float64}
    #no rotation escape hatch
    a ≈ b && return θ, ϕ
    #convert pole coordinates to unit cartesian points
    c₁ = sph2cart(a)
    c₂ = sph2cart(b)
    #angle between pole coordinates
    d = c₁ ⋅ c₂
    ψ = acos(c₁ ⋅ c₂)
    #axis of rotation
    k = unit(c₁ × c₂)
    #cartesian location of target point
    v = sph2cart(θ, ϕ)
    #rotate
    w = v*cos(ψ) + (k × v)*sin(ψ) + k*(k ⋅ v)*(1.0 - d)
    #convert back to spherical coordinates
    cart2usph(w...)
end

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
function agescaling(t)
    #expression for 1 Ga
    S₁ = 3.79e-14*(exp(6.93) - 1) + 5.84e-4
    #expression for t Ga
    Sₜ = 3.79e-14*(exp(6.93*t) - 1) + 5.84e-4*t
    #ratio
    Sₜ/S₁
end

function craterdensities(t)
    #mean crater radius for each bin [meters]
    r = 1e3*exp2.(𝐂["i"]/2 .+ 1/4)/2
    #frequency/density [craters/m^2]
    ρ = agescaling(t)*𝐂["N"]/1e6
    return r, ρ
end

function cratercounts(t, A)
    #radius bins and frequencies
    r, ρ = craterdensities(t)
    #counts [craters]
    n = ρ*A
    return r, ρ, n
end

#--------------------------------------
export Crater

struct Crater
    θ::Float64
    ϕ::Float64
    r::Float64
end

function Base.show(io::IO, c::Crater)
    θ, ϕ, r = map(x->round(x, sigdigits=4), (c.θ, c.ϕ, c.r))
    print(io, "crater θ=$θ, ϕ=$ϕ, r=$r")
end

#creates a randomly located crater with radius r
function Crater(r::Real)
    θ, ϕ = sphrand()
    Crater(θ, ϕ, r)
end

#creates a randomly located crater with radius r, using a specific random number generator
function Crater(r::Real, rng::AbstractRNG)
    θ, ϕ = sphrand(rng)
    Crater(θ, ϕ, r)
end

#multiplcation scales crater radius by a factor of f, returning a new Crater
function *(c::Crater, f::Real)
    @assert f >= 0 "crater radius cannot be negative"
    Crater(c.θ, c.ϕ, f*c.r)
end

#computes spherical distance between crater center and coordinate [θ,ϕ] using sphere radius r
sphdist(c::Crater, θ, ϕ, R=♂ᵣ) = sphdist(c.θ, c.ϕ, θ, ϕ, R)

#draws a circle on the sphere representing crater boundary
sphcirc(c::Crater, R=♂ᵣ; N::Int=75) = sphcirc(c.θ, c.ϕ, c.r, R; N=N)

#--------------------------------------
export GlobalPopulation

#==
This is a very low memory representation of a global crater population.
Crater coordinates are not stored, but produced on the fly when iterating.
A GlobalPopulation behaves like a random crater generator, so each
will have it's own random number generator for seeding and
reproducibility. Individual craters can only be generated via iteration,
not by indexing or anything else
==#
struct GlobalPopulation
    bins::Int64 #number of radius bins
    counts::Vector{Int64} #crater count in each bin
    N::Int64 #total number of craters
    r::Vector{Float64} #mean radius of each bin
    rng::Xoshiro
end

function Base.show(io::IO, P::GlobalPopulation)
    r = map(r->round(r, sigdigits=8), P.r)
    counts = [format(c, commas=true) for c in P.counts] 
    pretty_table(io,
        Any[1:P.bins r counts],
        ["", "Mean Radius [m]", "Count"]
    )
    N = format(P.N, commas=true)
    print(io, "Total craters: $N\n")
    print(io, "Random number generator seed: $(P.rng.seed)")
end

function GlobalPopulation(r::Vector{Float64}, counts::Vector{Int64}, seed=1)
    @assert length(r) == length(counts)
    GlobalPopulation(length(r), counts, sum(counts), r, Xoshiro(seed))
end

function GlobalPopulation(t::Real; rmin::Real=0, nmax::Real=Inf, seed=1)
    @assert t > 0
    @assert rmin >= 0
    @assert nmax >= 0
    #mean crater radius [m] and count for each bin
    r, _, n = cratercounts(t, ♂ₐ)
    #ROUND DOWN, just to be conservative
    counts = n .|> floor .|> Int64
    #don't use bins with huge counts for now (or zero craters)
    idx = @. (counts <= nmax) & (counts > 0) & (r >= rmin)
    r, counts = r[idx], counts[idx]
    #store it all in, effectively, a random crater generator
    return GlobalPopulation(r, counts, seed)
end

Base.eltype(::GlobalPopulation) = Crater

Base.length(P::GlobalPopulation) = P.N

function Base.iterate(P::GlobalPopulation,
                      state::NTuple{2,Int64}=(1,1)
                      )::Union{Nothing, Tuple{Crater,NTuple{2,Int64}}}
    @unpack bins, counts, r, rng = P
    #current bin and index within that bin
    bin, idx = state
    #stop at the end of the last bin
    bin > bins && return nothing
    #radius is always the same within a bin for this new crater
    @inbounds c = Crater(r[bin], rng)
    #continue in the same bin or start the next one
    @inbounds (idx > counts[bin]) ? (c, (bin+1,1)) : (c, (bin,idx+1))
end

#==============================================================================
The following functions and definitions handle the simulation of craters
impacting a hypothetical shoreline.
==============================================================================#

export SimulationResult
export segmentlengths

struct SimulationResult{T}
    #number of registered impacts
    impacts::Int64
    #fraction of shoreline that survived
    survived::Float64
    #fraction of shoreline destroyed
    destroyed::Float64
    #surviving shoreline segments
    segments::Vector{T}
    #all craters registered as impacting the line
    impactors::Vector{Crater}
end

function Base.show(io::IO, res::SimulationResult{T}) where {T}
    println(io, "SimulationResult{$T}")
    println(io, "  $(res.impacts) unique impacts")
    f = round(100*res.survived, sigdigits=6)
    println(io, "  $f % survived")
    f = round(100*res.destroyed, sigdigits=6)
    print(io, "  $f % destroyed")
end

function connected(a::NTuple{2,Float64}, b::NTuple{2,Float64})::Bool
    (@inbounds (a[1] == b[2]) | (b[1] == a[2])) ? true : false
end

# need to do this with a graph???
function segmentlengths(res::Vector{SphericalSegment},
                        θₛ::Float64, #segment latitude
                        R::Float64=♂ᵣ #sphere radius
                        )
    #pull out segments
    S = res.segments
    N = length(S)
    #fill an adjacency matrix
    A = spzeros(Bool, N, N)
    @inbounds for i ∈ 1:N, j ∈ 1:N
        if connected(S[i], S[j])
            A[i,j] = true
        end
    end
    return A
end

#computes segment lengths (in meters) of an impacted shoreline
function segmentlengths(segments::Vector{NTuple{2,Float64}},
                        θₛ::Float64, #segment latitude
                        R::Float64=♂ᵣ #sphere radius
                        )::Vector{Float64}

    #returned segment lengths
    seglen = Float64[]

    if length(S) == 1
        #a single segment should be a complete circle
        @assert S[1] == (0.0,𝛕)
        push!(seglen, 𝛕*R)
    else
        #multiple segments present
        for i ∈ 1:length(S)-1
            #segment length in meters
            push!(seglen, R*(S[i][2] - S[i][1]))
        end
        #final segment in radians
        Δϕ = S[end][2] - S[end][1]
        #check if it is distinct or wraps into the first seg
        if (S[1][1] == 0) & (S[end][2] == 𝛕)
            seglen[1] += R*Δϕ
        else
            push!(seglen, R*Δϕ)
        end
    end
    #segment distances need to be scaled by latitute
    seglen .*= sin(θₛ)

    return seglen
end

function segmentlengths(res::SimulationResult, θₛ::Real, R::Real=♂ᵣ)
    segmentlengths(res.segments, Float64(θₛ), Float64(R))
end

#--------------------------------------
#iso-latitude representative shoreline

export simulateimpacts

ℱᵣ(c::Crater, θ::Float64, Δϕ::Float64)::Float64 = sphdist(c, θ, c.ϕ + Δϕ) - c.r

function root(crater::Crater,
              θ::Float64,
              Δϕ₁::Float64,
              Δϕ₂::Float64,
              maxiter::Int64=1000)::Float64
    d₁ = ℱᵣ(crater, θ, Δϕ₁)
    d₂ = ℱᵣ(crater, θ, Δϕ₂)
    Δϕ = Inf
    δϕ = Inf
    d = Inf
    n::Int64 = 0
    #secant method with stringent termination tolerance
    while (abs(d₂ - d₁) > 1e-10) & (abs(δϕ) > 1e-10)
        #approximate root
        δϕ = d₁*(Δϕ₂ - Δϕ₁)/(d₂ - d₁)
        Δϕ = Δϕ₁ - δϕ
        d = ℱᵣ(crater, θ, Δϕ)
        #swaps
        Δϕ₁ = Δϕ₂
        Δϕ₂ = Δϕ
        d₁ = d₂
        d₂ = d
        #break on non-convergence
        n += 1
        n == maxiter && error("$maxiter iterations encoutered, Δϕ₁=$Δϕ₁, Δϕ₂=$Δϕ₂, δϕ=$δϕ, d₁=$d₁, d₂=$d₂, d=$d, crater=$crater")
    end
    return Δϕ
end

function intersection(crater::Crater, θₛ::Real, R::Real=♂ᵣ)
    @assert 0 <= θₛ <= π
    #crater parameters
    @unpack θ, ϕ, r = crater
    #double check that the crater overlaps the colatitude ring
    @assert R*abs(θ - θₛ) < r
    #find the intersection numerically/iteratively
    Δϕ = root(crater, θₛ, 0.0, π/1.1)
    #create a longitude interval with values ∈ [0,2π]
    ϕ₁, ϕ₂ = wrapangle(ϕ - Δϕ), wrapangle(ϕ + Δϕ)
    return ϕ₁, ϕ₂
end

function bite!(S::Vector{NTuple{2,T}}, sₙ::T, eₙ::T)::Nothing where {T<:Real}
    #number of stored intervals
    L = length(S)
    #check them all for partial or total removal
    i = 1
    stop = false
    while (i <= L) & !stop
        #stored interval
        s, e = S[i]
        #check various overlap cases
        if (s < sₙ) & (eₙ < e)
            #new interval is inside the stored one, split
            S[i] = (s, sₙ)
            insert!(S, i+1, (eₙ, e))
            stop = true #search is over b/c stored intervals don't overlap
        elseif (sₙ <= s) & (e <= eₙ)
            #new interval contains (or is identical to) the stored one, delete
            deleteat!(S, i)
            i -= 1
            L -= 1
        elseif (sₙ <= s) & (s < eₙ < e)
            #overlap on the lower side, crop up
            S[i] = (eₙ, e)
        elseif (s < sₙ < e) & (e <= eₙ)
            #overlap on the upper side, crop down
            S[i] = (s, sₙ)
        end
        i += 1
    end
end

function simulateimpacts(population::GlobalPopulation,
                         θₛ::Float64, #shoreline colatitude [0,π]
                         rₑ::Float64=1.0, #ejecta scaling of radius
                         Δ::Float64=0.0 #required overlap distance for impact to register
                         )::SimulationResult
    #check coordinate boundaries
    @assert 0.0 <= θₛ <= π "shoreline colatitude (θₛ) must be ∈ [0,π]"
    #check overlap distance
    @assert Δ >= 0.0 "overlap distance (Δ) must be positive"
    #start a shoreline to take bites out of
    segments = NTuple{2,Float64}[(0.0,𝛕)]
    #keep a list of craters that impact
    impactors = Crater[]
    #keep track of the total number of impacts
    n::Int64 = 0
    #now go through each crater, chopping up the shoreline as necessary
    for crater ∈ population
        #adjust radius for ejecta
        crater *= rₑ
        #short parameter names
        @unpack θ, ϕ, r = crater
        #distance from crater center to line
        dₛ = ♂ᵣ*abs(θₛ - θ)
        #check if the crater overlaps the line enough
        if dₛ < r - Δ
            #register the impact
            push!(impactors, crater)
            n += 1
            #find the longitude intersection interval
            ϕ₁, ϕ₂ = intersection(crater, θₛ)
            if (ϕ₂ < ϕ₁)# & (abs(ϕ₁ - ϕ₂) > π/6)
                bite!(segments, 0., min(ϕ₁, ϕ₂))
                bite!(segments, max(ϕ₁, ϕ₂), 𝛕)
            else
                bite!(segments, ϕ₁, ϕ₂)
            end
        end
    end
    #compute the fraction surviving
    f = sum(seg -> seg[2] - seg[1], segments)/𝛕
    #construct the whole result
    return SimulationResult(n, f, 1 - f, segments, impactors)
end

function simulateimpacts(t::Real, #time [Ga]
                         θₛ::Real, #colatitude of synthetic shoreline [rad]
                         rₑ::Real=1.0, #ejecta scaling of radius
                         Δ::Real=0.0; #required overlap distance for impact to register
                         rmin::Real=1e3, #smallest allowed crater radius [m]
                         nmax::Real=1_000_000, #maximum craters in bins, default small value
                         seed=1,
                         show::Bool=false)::SimulationResult
    #start up the crater population (impossible to have impacts where r < Δ)
    population = GlobalPopulation(t, rmin=max(rmin,Δ), nmax=nmax, seed=seed)
    #print the crater population table if desired
    show && println(population)
    #send the craters!
    simulateimpacts(population, Float64(θₛ), Float64(rₑ), Float64(Δ))
end

#--------------------------------------
#arbitrary segments

function simulateimpacts(population::GlobalPopulation,
                         segments::Vector{SphericalSegment},
                         rₑ::Float64=1.0,
                         Δ::Float64=0.0
                         )::SimulationResult
    #

end

end
