module ShorelineSurvival

import Base.*
using LinearAlgebra: ⋅, ×
using PrettyTables
using UnPack
using Random: AbstractRNG, MersenneTwister, seed!
using StaticArrays
using ForwardDiff: derivative
using Roots
using MultiAssign
using Formatting
using CSV
using DataFrames

const 𝛕 = 2π
export 𝛕

#==============================================================================
This section contains functions for doing various things in spherical geometry.
I use a colatitude coordinate θ (theta) ∈ [0,π] and a longitude
coordinate ϕ (phi) ∈ [0,2π].
==============================================================================#

#--------------------------------------
export sphrand

function sphrand(rng::AbstractRNG)::NTuple{2,Float64}
    θ = acos(1 - 2*rand(rng))
    ϕ = 𝛕*rand(rng)
    return θ, ϕ
end

function sphrand()::NTuple{2,Float64}
    θ = acos(1 - 2*rand())
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
export latlon2sph

function latlon2sph(lat::Real, lon::Real)
    @assert -90 <= lat <= 90 "latitude must be ∈ [-90,90]"
    @assert -180 <= lat <= 180 "longitude must be ∈ [-180,180]"
    θ = -lat*(π/180) + π/2
    ϕ = lon*(π/180) + π
    return θ, ϕ
end

function latlon2sph(lat::AbstractVector{T}, lon::AbstractVector{T}) where {T}
    @assert length(lat) == length(lon)
    @multiassign θ, ϕ = zeros(T, length(lat))
    @inbounds for i ∈ 1:length(lat)
        θ[i], ϕ[i] = latlon2sph(lat[i], lon[i])
    end
    return θ, ϕ
end

#--------------------------------------
export sph2cart

#assumes radius is 1
function sph2cart(θ::T, ϕ::T) where {T}
    sₜ, cₜ = sincos(θ)
    sₚ, cₚ = sincos(ϕ)
    return SVector{3,T}(sₜ*cₚ, sₜ*sₚ, cₜ)
end

sph2cart(θ, ϕ, r) = r*sph2cart(θ, ϕ)

function sph2cart(θ::AbstractVector{T}, ϕ::AbstractVector{T}, r::T) where {T}
    @assert length(θ) == length(ϕ)
    @multiassign x, y, z = similar(θ)
    @inbounds for i ∈ 1:length(x)
        x[i], y[i], z[i] = sph2cart(θ[i], ϕ[i], r)
    end
    return x, y, z
end

function sph2cart(θ::AbstractVector{T}, ϕ::AbstractVector{T}) where {T}
    sph2cart(θ, ϕ, one(T))
end

#--------------------------------------
export cart2sph, cart2usph

function cart2sph(x::T, y::T, z::T) where {T}
    r = sqrt(x*x + y*y + z*z)
    θ = acos(z/r)
    ϕ = ↻(atan(y,x))
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
function cart2usph(x::T, y::T, z::T) where {T}
    θ, ϕ, _ = cart2sph(x, y, z)
    return θ, ϕ
end

cart2usph(v::SVector{3,T}) where {T} = cart2usph(v...)

#--------------------------------------
#arc lengths and spherical distances

export arclength, sphdist

#assumes vectors have length 1
function arclength(c₁::SVector{3,T}, c₂::SVector{3,T}) where {T}
    (c₁ == c₂) && return zero(T)
    acos(c₁ ⋅ c₂)
end

function arclength(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where {T}
    acos(sph2cart(θ₁, ϕ₁) ⋅ sph2cart(θ₂, ϕ₂))
end

sphdist(θ₁, ϕ₁, θ₂, ϕ₂, R) = R*arclength(θ₁, ϕ₁, θ₂, ϕ₂)

#--------------------------------------
#misc

export ↻, unit, sphcirc, unitnormal

#wraps an angle into [0,2π] and appears to be quicker than using remainder
function ↻(θ)
    while θ < 0; θ += 𝛕; end
    while θ > 𝛕; θ -= 𝛕; end
    return θ
end

function unit(v::SVector{3,T}) where {T}
    x, y, z = v
    L = sqrt(x*x + y*y + z*z)
    return SVector{3,T}(x/L, y/L, z/L)
end

unitnormal(a::SVector{3,T}, b::SVector{3,T}) where {T} = unit(a × b)

function sphcirc(θ::T, ϕ::T, r::T, R=♂ᵣ; N::Int=50) where {T}
    #vector from center of sphere to center of circle
    C = sph2cart(θ, ϕ, convert(T, R))
    #unit vector from sphere center to circle center, normal to circle's plane
    n = unit(C)
    #unit vector perpendicular to n in the x,y plane
    u = SVector{3,T}(-sin(ϕ), cos(ϕ), 0.0)
    #unit vector perpendicular to both n and u using cross product
    v = n × u
    #crater radius arc length
    𝓁 = r/R
    #reduced cartesian distance from origin
    D = C*cos(𝓁)
    #reduced cartesian radius of circle for curvature
    d = R*sin(𝓁)
    #create vectors of coordinates representing the circle
    @multiassign x, y, z = zeros(T, N)
    u₁, u₂, u₃ = u
    v₁, v₂, v₃ = v
    D₁, D₂, D₃ = D
    @inbounds for (i,ψ) ∈ enumerate(LinRange(0, 𝛕, N))
        s = sin(ψ)
        c = cos(ψ)
        x[i] = D₁ + d*(s*u₁ + c*v₁)
        y[i] = D₂ + d*(s*u₂ + c*v₂)
        z[i] = D₃ + d*(s*u₃ + c*v₃)
    end
    return x, y, z
end

function checkcoord(θ, ϕ)::Nothing
    @assert 0.0 <= θ <= π
    @assert 0.0 <= ϕ <= 𝛕
    nothing
end

#==============================================================================
Here two simple types for spherical geometry are defined, along with
some basic operations on them as wrappers of the functions above.
==============================================================================#

#--------------------------------------
export SphericalPoint

struct SphericalPoint{T}
    θ::T
    ϕ::T
end

#Base.show(io::IO, p::SphericalPoint) = print(io, "(θ=$(p.θ), ϕ=$(p.ϕ))")

SphericalPoint(x::NTuple{2,T}) where {T} = @inbounds SphericalPoint{T}(x[1], x[2])

sph2cart(p::SphericalPoint) = sph2cart(p.θ, p.ϕ)

arclength(a::SphericalPoint{T}, b::SphericalPoint{T}) where {T} = arclength(a.θ, a.ϕ, b.θ, b.ϕ)

checkpoint(p::SphericalPoint)::Nothing = checkcoord(p.θ, p.ϕ)

#--------------------------------------
export SphericalSegment
export minplanecolat

struct SphericalSegment{T}
    a::SphericalPoint{T}
    b::SphericalPoint{T}
end

function SphericalSegment(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where {T}
    SphericalSegment{T}(
        SphericalPoint(θ₁, ϕ₁),
        SphericalPoint(θ₂, ϕ₂)
    )
end

function SphericalSegment(a::NTuple{2,T}, b::NTuple{2,T}) where {T}
    SphericalSegment{T}(SphericalPoint(a), SphericalPoint(b))
end

arclength(s::SphericalSegment) = arclength(s.a, s.b)

sphdist(s::SphericalSegment, R::Real=♂ᵣ) = R*arclength(s)

sph2cart(s::SphericalSegment) = sph2cart(s.a), sph2cart(s.b)

unitnormal(s::SphericalSegment) = unitnormal(sph2cart(s.a), sph2cart(s.b))

function checksegment(s::SphericalSegment, maxarc=π/6)::Nothing
    checkpoint(s.a)
    checkpoint(s.b)
    𝓁 = arclength(s)
    if 𝓁 > maxarc
        p = round(100*𝓁/𝛕, sigdigits=4)
        error("unusually large segment with arclength=$𝓁 or ~$p % of 2π")
    end
    nothing
end

function commonendpoint(s₁::SphericalSegment{T}, s₂::SphericalSegment{T})::Bool where {T}
    c₁ = (sph2cart(s₁.a), sph2cart(s₁.b))
    c₂ = (sph2cart(s₂.a), sph2cart(s₂.b))
    for p₁ ∈ c₁, p₂ ∈ c₂
        (p₁ == p₂) && return true
    end
    return false
end

#--------------------------------------
export CartesianSegment

struct CartesianSegment{T}
    a::SVector{3,T}
    b::SVector{3,T}
end

function CartesianSegment(s::SphericalSegment{T}) where {T}
    CartesianSegment{T}(sph2cart(s.a), sph2cart(s.b))
end

function minplanecolat(c::CartesianSegment)
    n = unitnormal(c.a, c.b)
    z = @inbounds abs(n[3])
    return asin(z)
end

#==============================================================================
This stuff is for rotating points/segments over the sphere
==============================================================================#

export SphericalRotation, identityrotation, rotate, unrotate

struct SphericalRotation{T}
    d::T
    ψ::T
    k::SVector{3,T}
end

function identityrotation(T::Type=Float64)
    o, z = one(T), zero(T)
    SphericalRotation(o, z, unit(SVector{3,T}(o, o, o)))
end

function SphericalRotation(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where {T}
    #coordinate integrity
    checkcoord(θ₁, ϕ₁)
    checkcoord(θ₂, ϕ₂)
    #cartesian vectors
    c₁ = sph2cart(θ₁, ϕ₁)
    c₂ = sph2cart(θ₂, ϕ₂)
    #account for unlikely cases with zero rotation to prevent NaN
    (c₁ == c₂) && return identityrotation(T)
    #angle between pole coordinates
    d = c₁ ⋅ c₂
    ψ = acos(c₁ ⋅ c₂)
    #axis of rotation
    k = unit(c₁ × c₂)
    #construct
    SphericalRotation(d, ψ, k)
end

#--------------------------------------

function rotate(v::SVector{3,T}, d::T, ψ::T, k::SVector{3,T}) where {T}
    s, c = sincos(ψ)
    v*c + (k × v)*s + k*(k ⋅ v)*(1 - d)
end

function rotate(v::SVector{3,T}, ρ::SphericalRotation{T}) where {T}
    rotate(v, ρ.d, ρ.ψ, ρ.k)
end

function rotate(c::CartesianSegment{T}, ρ::SphericalRotation{T}) where {T}
    CartesianSegment{T}(rotate(c.a, ρ), rotate(c.b, ρ))
end

function rotate(θ::T, ϕ::T, ρ::SphericalRotation{T}) where {T}
    cart2usph(rotate(sph2cart(θ, ϕ), ρ))
end

function rotate(p::SphericalPoint{T}, ρ::SphericalRotation{T}) where {T}
    SphericalPoint{T}(rotate(p.θ, p.ϕ, ρ))
end

function rotate(s::SphericalSegment{T}, ρ::SphericalRotation{T}) where {T}
    SphericalSegment{T}(rotate(s.a, ρ), rotate(s.b, ρ))
end

#--------------------------------------

function unrotate(θ::T, ϕ::T, ρ::SphericalRotation{T}) where {T}
    #to reverse the rotation, k changes sign
    cart2usph(rotate(sph2cart(θ, ϕ), ρ.d, ρ.ψ, -ρ.k))
end

function unrotate(p::SphericalPoint{T}, ρ::SphericalRotation{T}) where {T}
    SphericalPoint(unrotate(p.θ, p.ϕ, ρ))
end

function unrotate(s::SphericalSegment{T}, ρ::SphericalRotation{T}) where {T}
    SphericalSegment(unrotate(s.a, ρ), unrotate(s.b, ρ))
end

#==============================================================================
This type sets up a parameterized equation for a great circle through two
points, with periodic parameter range ∈ [0,2π]
# https://math.stackexchange.com/questions/1783746/equation-of-a-great-circle-passing-through-two-points
==============================================================================#

export GreatCircle
export sph, colat
export 𝒟, 𝒟′, 𝒟′′, argmindist, argmincolat

struct GreatCircle{T}
    v::SVector{3,T}
    w::SVector{3,T}
    function GreatCircle(v₁::SVector{3,T}, v₂::SVector{3,T}) where {T}
        @assert !isapprox(v₁, v₂, rtol=1e-12)
        d = v₁ ⋅ v₂
        f = sqrt(1 - d^2)
        α = -d/f
        β = 1/f
        w = α*v₁ + β*v₂
        new{T}(v₁, w)
    end
end

GreatCircle(c::CartesianSegment) = GreatCircle(c.a, c.b)

function GreatCircle(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where {T}
    GreatCircle(sph2cart(θ₁, ϕ₁), sph2cart(θ₂, ϕ₂))
end

function GreatCircle(a::SphericalPoint{T}, b::SphericalPoint{T}) where {T}
    GreatCircle(a.θ, a.ϕ, b.θ, b.ϕ)
end

GreatCircle(s::SphericalSegment) = GreatCircle(s.a, s.b)

#t is a parameter ∈ [0,2π] defining the circle
function (C::GreatCircle)(t)
    s, c = sincos(t)
    c*C.v + s*C.w
end

sph(C::GreatCircle, t) = SphericalPoint(cart2usph(C(t)))

function colat(C::GreatCircle, t)
    x, y, z = C(t)
    θ, _ = cart2sph(x, y, z)
    return θ
end

#objective function for distance between point and great circle
𝒟(C::GreatCircle, t, X) = sum((C(t) - X).^2)

𝒟′(C::GreatCircle, t, v) = derivative(x->𝒟(C,x,v), t)

𝒟′′(C::GreatCircle, t, v) = derivative(x->𝒟′(C,x,v), t)

function argmindist(C::GreatCircle{T},
                    p::SphericalPoint{T},
                    tol::Float64=1e-14,
                    maxiter::Int64=1000) where {T}
    #cartesian location of target point
    v = sph2cart(p)
    #newton's minimization
    t = abs(𝒟′(C, 0.0, v)) < abs(𝒟′(C, 1.0, v)) ? 0.0 : 1.0
    Δ = Inf
    Δₚ = 0.0
    f′ = 𝒟′(C, t, v)
    f′′ = 𝒟′′(C, t, v)
    n::Int64 = 0
    while (abs(Δ - Δₚ) > tol) | (abs(f′) > tol)
        Δ, Δₚ = f′/f′′, Δ
        t = ↻(t - Δ)
        f′ = 𝒟′(C, t, v)
        f′′ = 𝒟′′(C, t, v)
        n += 1
        (n == maxiter) && error("$maxiter iterations encountered $Δ $Δₚ $f′ $f′′ $C")
    end
    𝒟(C, t, v) < 𝒟(C, t + π, v) ? t : ↻(t + π)
end

argmincolat(C::GreatCircle{T}) where {T} = argmindist(C, SphericalPoint(zero(T), zero(T)))

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

struct Crater{T}
    θ::T
    ϕ::T
    r::T
end

#function Base.show(io::IO, c::Crater)
#    θ, ϕ, r = map(x->round(x, sigdigits=4), (c.θ, c.ϕ, c.r))
#    print(io, "crater θ=$θ, ϕ=$ϕ, r=$r")
#end

#creates a randomly located crater with radius r
function Crater(r::T) where {T}
    θ, ϕ = sphrand()
    Crater(T(θ), T(ϕ), r)
end

#creates a randomly located crater with radius r, using a specific random number generator
function Crater(r::T, rng::AbstractRNG) where {T}
    θ, ϕ = sphrand(rng)
    Crater(T(θ), T(ϕ), r)
end

#multiplcation scales crater radius by a factor of f, returning a new Crater
function *(c::Crater{T}, f::Real) where {T}
    @assert f >= 0 "crater radius cannot be negative"
    Crater(c.θ, c.ϕ, T(f*c.r))
end

#computes spherical distance between crater center and coordinate [θ,ϕ] using sphere radius r
sphdist(c::Crater, θ, ϕ, R=♂ᵣ) = sphdist(c.θ, c.ϕ, θ, ϕ, R)

#draws a circle on the sphere representing crater boundary
sphcirc(c::Crater, R=♂ᵣ; N::Int=50) = sphcirc(c.θ, c.ϕ, c.r, R; N=N)

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
struct GlobalPopulation
    bins::Int64 #number of radius bins
    counts::Vector{Int64} #crater count in each bin
    N::Int64 #total number of craters
    r::Vector{Float64} #mean radius of each bin
    rng::MersenneTwister
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
end

function GlobalPopulation(r::Vector, counts::Vector{Int64}, seed=1) where {T}
    @assert length(r) == length(counts)
    GlobalPopulation(length(r), counts, sum(counts), r, MersenneTwister(seed))
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
    #iterate more efficiently with larger craters in the first bins
    idx = sortperm(r, rev=true)
    r = r[idx]
    counts = counts[idx]
    #store it all in, effectively, a random crater generator
    return GlobalPopulation(r, counts, seed)
end

Base.eltype(::GlobalPopulation) = Crater{Float64}

Base.length(P::GlobalPopulation) = P.N

function Base.iterate(P::GlobalPopulation,
                      state::NTuple{2,Int64}=(1,1)
                      )::Union{Nothing,Tuple{Crater{Float64},NTuple{2,Int64}}}
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
export segmentdistances

struct SimulationResult{T}
    impacts::Int64 #number of registered impacts
    survived::Float64 #fraction of shoreline that survived
    destroyed::Float64 #fraction of shoreline destroyed
    segments::Vector{T} #surviving shoreline segments
    impactors::Vector{Crater} #all craters registered as impacting the line
end

function Base.show(io::IO, res::SimulationResult{T}) where {T}
    println(io, "SimulationResult{$T}")
    println(io, "  $(res.impacts) unique impacts")
    f = round(100*res.survived, sigdigits=6)
    println(io, "  $f % survived")
    f = round(100*res.destroyed, sigdigits=6)
    print(io, "  $f % destroyed")
end

#computes segment lengths (in meters) of an impacted shoreline
function segmentdistances(S::Vector{NTuple{2,Float64}},
                          θₛ::Float64, #segment latitude
                          R::Float64=♂ᵣ #sphere radius
                          )::Vector{Float64}
    if length(S) == 1
        a = [𝛕]
    else
        a = map(s->s[2]-s[1], S)
        #check if first and last segments actually wrap
        if (S[1][1] == 0) & (S[end][2] == 𝛕)
            a[1] += pop!(a)
        end
    end
    #scale by radius and latitude to get distance in meters
    return a*R*sin(θₛ)
end

function segmentdistances(S::Vector{SphericalSegment{T}}, R::Float64=♂ᵣ) where {T}
    #assume the segments are in order
    a = arclength.(S)
    𝓁 = T[a[1]]
    for i ∈ 2:length(S)
        if commonendpoint(S[i], S[i-1])
            𝓁[end] += a[i]
        else
            push!(𝓁, a[i])
        end
    end
    #handle possible wrapping
    if commonendpoint(S[1], S[end])
        𝓁[1] += pop!(𝓁)
    end
    #remember to apply the radius
    return R*𝓁
end

function segmentdistances(res::SimulationResult, args...)
    segmentdistances(res.segments, args...)
end

#--------------------------------------
#iso-latitude representative shoreline

function lonshift(crater::Crater,
                  θ::Float64,
                  Δϕ₁::Float64,
                  Δϕ₂::Float64,
                  tol::Float64=1e-12,
                  maxiter::Int64=1000)::Float64
    #cartesian crater center
    χ = sph2cart(crater.θ, crater.ϕ)
    #arclength of crater radius
    𝓁₀ = crater.r/♂ᵣ
    #find where the crater meets θ 
    find_zero(
        Δϕ->arclength(sph2cart(θ, crater.ϕ + Δϕ), χ) - 𝓁₀,
        (Δϕ₁, Δϕ₂),
        Roots.Order0(),
        atol=tol,
        rtol=0.0,
        xatol=tol,
        xrtol=0.0,
        maxevals=maxiter,
        strict=true
    )
end

function intersection(crater::Crater{Float64}, θₛ::Real, R::Real=♂ᵣ)
    #crater parameters
    @unpack θ, ϕ, r = crater
    #double check that the crater overlaps the colatitude ring
    @assert R*abs(θ - θₛ) < r
    #find the intersection numerically/iteratively
    Δϕ = lonshift(crater, θₛ, 0.0, π/3)
    #create a longitude interval with values ∈ [0,2π]
    ϕ₁, ϕ₂ = ↻(ϕ - Δϕ), ↻(ϕ + Δϕ)
    return ϕ₁, ϕ₂
end

function overlapcase(s::T, e::T, sₙ::T, eₙ::T)::Int64 where {T}
    if (sₙ >= e) | (eₙ <= s)
        #no overlap
        return 0
    elseif (s < sₙ) & (eₙ < e)
        #new interval is inside
        return 1
    elseif (sₙ <= s) & (e <= eₙ)
        #overlap on the lower side
        return 2
    elseif (sₙ <= s) & (s < eₙ < e)
        #overlap on the upper side
        return 3
    elseif (s < sₙ < e) & (e <= eₙ)
        return 4
    else
        error("overlap case failure")
    end
end

function clip!(S::Vector{NTuple{2,T}}, sₙ::T, eₙ::T)::Bool where {T}
    #number of stored intervals
    L = length(S)
    #check them all for partial or total removal
    i = 1
    stop = false
    impacted = false
    while (i <= L) & !stop
        #stored interval
        s, e = S[i]
        #get overlap type (or non-overlap)
        case = overlapcase(s, e, sₙ, eₙ)
        #check various overlap cases
        if case == 1
            #new interval is inside the stored one, split
            S[i] = (s, sₙ)
            insert!(S, i+1, (eₙ, e))
            stop = true #search is over b/c stored intervals don't overlap
            impacted = true
        elseif case == 2
            #new interval contains (or is identical to) the stored one, delete
            deleteat!(S, i)
            i -= 1
            L -= 1
            impacted = true
        elseif case == 3
            #overlap on the lower side, crop up
            S[i] = (eₙ, e)
            impacted = true
        elseif case == 4
            #overlap on the upper side, crop down
            S[i] = (s, sₙ)
            impacted = true
        end
        i += 1
    end
    return impacted
end

function simulateimpacts(population::GlobalPopulation,
                         θₛ::Float64, #shoreline colatitude [0,π]
                         rₑ::Float64, #ejecta scaling of radius
                         Δ::Float64 #required overlap distance for impact to register
                         )::SimulationResult
    #check coordinate boundaries
    @assert 0.0 <= θₛ <= π "shoreline colatitude (θₛ) must be ∈ [0,π]"
    #start a shoreline to take bites out of
    segs = NTuple{2,Float64}[(0.0,𝛕)]
    #store craters that impact
    impactors = Set{Crater}()
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
            #find the longitude intersection interval
            ϕ₁, ϕ₂ = intersection(crater, θₛ)
            #clip overlapping portions
            if ϕ₂ < ϕ₁ #intersection interval wraps over 2π
                if clip!(segs, 0., min(ϕ₁, ϕ₂)) | clip!(segs, max(ϕ₁, ϕ₂), 𝛕)
                    push!(impactors, crater)
                end
            else
                clip!(segs, ϕ₁, ϕ₂) && push!(impactors, crater)
            end
        end
    end
    #compute the fraction surviving
    f = (length(segs) > 1) ? sum(x->x[2]-x[1], segs)/𝛕 : 0.0
    #construct the final result
    SimulationResult(
        length(impactors),
        f,
        1 - f,
        segs,
        collect(impactors)
    )
end

function simulateimpacts(population::GlobalPopulation,
                         θₛ::Real,
                         rₑ::Real,
                         Δ::Real)::SimulationResult
    simulateimpacts(population, Float64(θₛ), Float64(rₑ), Float64(Δ))
end

#--------------------------------------
#arbitrary segment simulations for mapped shorelines

export readsegments

#assumes lat ∈ [-90,90] and lon ∈ [-180,180]
function readsegments(fn::String;
                      minarc::Real=0.0,
                      lonname::String="lon",
                      latname::String="lat",
                      T::Type=Float64)
    @assert minarc >= 0.0
    #read the table
    df = CSV.read(fn, DataFrame)
    lat = df[!,latname]
    lon = df[!,lonname]
    #convert to radians
    θ, ϕ = latlon2sph(collect(T,lat), collect(T,lon))
    N = length(θ)
    L = N - 1
    #accumulate segments
    S = SphericalSegment{T}[]
    i = 1
    while i <= L
        #accumulate distance until exceeding mindist
        d = zero(T)
        j = i
        while (d <= minarc) & (j < N)
            d = arclength(θ[i], ϕ[i], θ[j+1], ϕ[j+1])
            j += 1
        end
        #add a segment
        push!(S, SphericalSegment(θ[i], ϕ[i], θ[j], ϕ[j]))
        #set the new starting point
        i = j
    end
    return S
end

function colatrange(S::Vector{SphericalSegment{T}}) where {T}
    θa = map(s->s.a.θ, S)
    θb = map(s->s.b.θ, S)
    θmin = min(minimum(θa), minimum(θb))
    θmax = max(maximum(θa), maximum(θb))
    return θmin, θmax
end

function intersection(C::GreatCircle{T},
                      θ::T, #target colatitude
                      tₘ::T, #parameter of minimum colatitude
                      t₀::T, #other starting point
                      maxiter::Int64=1000,
                      tol::T=1e-12) where {T}
    find_zero(
        x->colat(C, x) - θ,
        (tₘ, t₀),
        Roots.A42(),
        atol=tol,
        rtol=zero(T),
        xatol=tol,
        xrtol=zero(T),
        strict=true,
        maxevals=maxiter
    )
end

#parameter arguments where C intersects colatitude
function intersections(C::GreatCircle, θ::Float64)::NTuple{2,Float64}
    #minimize the colatitude along the great circle
    tₘ = argmincolat(C)
    #first point where great circle meets θ
    t₁ = intersection(C, θ, tₘ, tₘ - π/4)
    #the other point where they meet, by symmetry
    t₂ = tₘ + (tₘ - t₁)
    #t₂ may be > 2π, but this is handled later
    return t₁, t₂
end

function newseg(C::GreatCircle,
                t₁::Float64,
                t₂::Float64,
                ρ::SphericalRotation)::SphericalSegment
    unrotate(SphericalSegment(sph(C, t₁), sph(C, t₂)), ρ)
end

function clip!(segs,
               csegs,
               osegs,
               i::Int64,
               C::GreatCircle{T},
               ρ::SphericalRotation{T},
               𝓁ᵢ::Float64,
               sₙ::Float64,
               eₙ::Float64)::Tuple{Int64,Bool} where {T}
    #==========================================================
    If the intersection interval overruns 2π, the whole thing
    gets wrapped to the beginning. This is ok because the
    segment interval always starts at 0.0 and should never be
    very large (as check by s) and craters should not be
    so large that intersections span hemispheres or anything
    crazy like that.
    ==========================================================#
    if eₙ > 𝛕
        sₙ -= 𝛕
        eₙ -= 𝛕
    end
    #always have the parameter arguments of the segment's end points
    s, e = 0.0, 𝓁ᵢ
    #check various overlap cases
    impacted = true
    ΔL = 0
    case = overlapcase(s, e, sₙ, eₙ)
    if case == 1
        #split the segment, first part
        n = newseg(C, s, sₙ, ρ)
        segs[i] = newseg(C, s, sₙ, ρ)
        csegs[i] = CartesianSegment(segs[i])
        #second part
        n = newseg(C, eₙ, e, ρ)
        insert!(segs, i+1, n)
        insert!(csegs, i+1, CartesianSegment(segs[i+1]))
        insert!(osegs, i+1, osegs[i])
        ΔL += 1
    elseif case == 2
        #intersection contains the segment, discard the seg
        deleteat!(segs, i)
        deleteat!(csegs, i)
        deleteat!(osegs, i)
        ΔL = -1
    elseif case == 3
        #overlap on the lower side, crop up
        n = newseg(C, eₙ, e, ρ)
        segs[i] = n
        csegs[i] = CartesianSegment(segs[i])
    elseif case == 4
        #overlap on the upper side, crop down
        n = newseg(C, s, sₙ, ρ)
        segs[i] = n
        csegs[i] = CartesianSegment(segs[i])
    else
        impacted = false
    end
    return ΔL, impacted
end

function simulateimpacts(population::GlobalPopulation,
                         segs::Vector{SphericalSegment{𝒯}},
                         rₑ::Float64=1.0,
                         Δ::Float64=0.0,
                         minarc::Float64=1/♂ᵣ) where {𝒯}
    #check over segment coordinates
    foreach(checksegment, segs)
    L = length(segs)
    #find latitude range of segments
    θmin, θmax = colatrange(segs)
    #store initial sum of segment arclengths to compare with
    A₀ = sum(map(arclength, segs))
    #make a copy of the segments before taking bites out of them
    ssegs::Vector{SphericalSegment{𝒯}} = deepcopy(segs)
    #keep a cartesian mirror of the segments to speed up first filter
    csegs::Vector{CartesianSegment{𝒯}} = map(CartesianSegment, segs)
    #finally, copy of original segments for most accurate great circles
    osegs::Vector{CartesianSegment{𝒯}} = deepcopy(csegs)
    #store craters that impact
    impactors = Set{Crater}()
    #the arclength of the overlap buffer
    Δᵣ = Δ/♂ᵣ
    #iterate through the entire crater population
    for (count, crater) ∈ enumerate(population)
        #adjust crater radius for ejecta and unpack
        @unpack θ, ϕ, r = crater*rₑ
        #cartesian crater center
        χ = sph2cart(θ, ϕ) 
        #arclength of crater radius
        𝓁ᵣ = r/♂ᵣ
        #============================================================
        The first check for intersection is simply whether the crater
        is so far from the colatitude range of the putative shoreline
        segments that it's impossible for it to touch any of them
        ============================================================#
        if (θmin - 𝓁ᵣ) <= θ <= (θmax + 𝓁ᵣ)
            #set up rotation of this crater to the north pole
            ρ = SphericalRotation(θ, ϕ, zero(𝒯), zero(𝒯))
            #every segment has to be checked sadly
            i = 1
            while i <= L
                @inbounds cᵢ = csegs[i]
                #========================================================
                An easy second check is whether the distance/arclength
                between crater center and segment endpoints far exceeds
                the crater's radius
                ========================================================#
                if (acos(χ ⋅ cᵢ.a) - 𝓁ᵣ < π/4) & (acos(χ ⋅ cᵢ.b) - 𝓁ᵣ < π/4)
                    #original segment from which cᵢ may be a portion
                    @inbounds oᵢ = osegs[i]
                    #rotate to put crater center at the north pole
                    v = rotate(oᵢ, ρ)
                    #====================================================
                    The third check is whether the plane formed
                    by the rotated segment runs through the crater,
                    which is now at the north pole. This is relatively
                    easy to check and should save a fair amount of time.
                    This is simultaneously a check that the overlap meets
                    the minimum requirement Δ.
                    ====================================================#
                    if minplanecolat(v) + Δᵣ <= 𝓁ᵣ
                        #================================================
                        By this stage optimization doesn't matter much 
                        because the bulk of the work is done rejecting
                        intersections before this branch is reached.
                        Things still need to be robust, of course.
                        ================================================#
                        #arclength of the actual segment
                        𝓁ᵢ = @inbounds arclength(ssegs[i])
                        #check if the segment is too small to count
                        if 𝓁ᵢ < minarc
                            #remove it
                            deleteat!(ssegs, i)
                            deleteat!(csegs, i)
                            deleteat!(osegs, i)
                            L -= 1
                            i -= 1
                        else
                            #paramaterize the rotated segment's great circle
                            C = GreatCircle(v)
                            #parameter values where C intersects the crater
                            t₁, t₂ = intersections(C, 𝓁ᵣ)
                            #sanity check that intersection segment is not larger than crater
                            Δt = abs(t₂ - t₁)
                            d = 2𝓁ᵣ
                            δ = (Δt - d)/d
                            if (Δt > d) & (δ > 1e-4)
                                println(stdout, arclength(segs[i]))
                                println(stdout, rotate(segs[i], rotation))
                                println(stdout, 𝓁ᵣ, '\n')
                                flush(stdout)
                                #error("intersection interval $Δt larger than crater $d by $δ")
                            end
                            #now check for genuine overlap
                            ΔL, impacted = clip!(ssegs, csegs, osegs, i, C, ρ, 𝓁ᵢ, t₁, t₂)
                            #we have impact!
                            if impacted
                                #register the crater
                                push!(impactors, crater)
                                #handle possible changes in list length
                                L += ΔL
                                i += ΔL
                            end
                        end
                    end
                end
                i += 1
            end
        end
        #occasionally refresh the colatitude range
        if count % 1000 == 0
            θmin, θmax = colatrange(ssegs)
        end
    end
    #final sum of segment arclengths
    A = sum(map(arclength, ssegs))
    #fraction surviving
    f = A/A₀
    #final constructionu
    SimulationResult(
        length(impactors),
        f,
        1 - f,
        ssegs,
        collect(impactors)
    )
end

function simulateimpacts(population::GlobalPopulation,
                         segments::Vector{SphericalSegment{T}},
                         rₑ::Real,
                         Δ::Real) where {T}
    simulateimpacts(population, segments, T(rₑ), T(Δ))
end

#--------------------------------------
#convenience and barrier function

export simulateimpacts

function simulateimpacts(t::Real, #time [Ga]
                         shoreline, #putative shoreline segments or latitude
                         rₑ::Real=1.0, #ejecta scaling of radius
                         Δ::Real=0.0; #required overlap distance for impact to register                         rmin::Real=1e3, #smallest allowed crater radius [m]
                         rmin::Real=1e3, #smallest allowed crater radius [m]
                         nmax::Real=1_000_000, #maximum craters in bins, default small value
                         seed=1,
                         show::Bool=false)::SimulationResult
    #check overlap distance
    @assert Δ >= 0 "overlap distance (Δ) must be positive"
    #check ejecta radius multiple
    @assert rₑ >= 1 "ejecta radius multiple (rₑ) must be greater than or equal to 1"
    #start up the crater population (impossible to have impacts where r < Δ)
    population = GlobalPopulation(t, rmin=max(rmin,Δ), nmax=nmax, seed=seed)
    #print the crater population table if desired
    show && println(population)
    #send the craters!
    simulateimpacts(population, shoreline, rₑ, Δ)
end

end
