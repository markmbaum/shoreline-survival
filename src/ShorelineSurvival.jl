module ShorelineSurvival

import Base.*
using LinearAlgebra: ⋅, ×
using PrettyTables
using UnPack
using Random: AbstractRNG, Xoshiro, seed!
using StaticArrays
using ForwardDiff: derivative
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

function latlon2sph(lat::AbstractVector{T}, lon::AbstractVector{T}) where {T<:Real}
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
function sph2cart(θ::T, ϕ::T) where {T<:Real}
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

function sph2cart(θ::AbstractVector{T}, ϕ::AbstractVector{T}) where {T<:Real}
    sph2cart(θ, ϕ, one(T))
end

#--------------------------------------
export cart2sph, cart2usph

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
export arclength, sphdist, unit, sphcirc, wrapangle, unitnormal

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

function unitnormal(a::SVector{3,T}, b::SVector{3,T}) where {T<:Real}
    unit(a × b)
end

function sphcirc(θ::T, ϕ::T, r::T, R=♂ᵣ; N::Int=50) where {T<:Real}
    #vector from center of sphere to center of circle
    C = sph2cart(θ, ϕ, convert(T, R))
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

struct SphericalPoint
    θ::Float64
    ϕ::Float64
end

Base.show(io::IO, p::SphericalPoint) = print(io, "(θ=$(p.θ), ϕ=$(p.ϕ))")

SphericalPoint(x::NTuple{2}) = @inbounds SphericalPoint(x[1], x[2])

function Base.isapprox(a::SphericalPoint, b::SphericalPoint)::Bool
    (a.θ ≈ b.θ) & (a.ϕ ≈ b.ϕ)
end

sph2cart(p::SphericalPoint)::SVector{3,Float64} = sph2cart(p.θ, p.ϕ)

function arclength(a::SphericalPoint, b::SphericalPoint)::Float64
    arclength(a.θ, a.ϕ, b.θ, b.ϕ)
end

checkpoint(p::SphericalPoint)::Nothing = checkcoord(p.θ, p.ϕ)

#--------------------------------------
export SphericalSegment

struct SphericalSegment
    a::SphericalPoint
    b::SphericalPoint
end

function Base.show(io::IO, s::SphericalSegment)
    print(io, "SphericalSegment:\n  a=$(s.a)\n  b=$(s.b)")
end

function SphericalSegment(θ₁, ϕ₁, θ₂, ϕ₂)
    SphericalSegment(
        SphericalPoint(θ₁, ϕ₁),
        SphericalPoint(θ₂, ϕ₂)
    )
end

function SphericalSegment(a::NTuple{2,Float64}, b::NTuple{2,Float64})
    @inbounds SphericalSegment(
        SphericalPoint(a[1], a[2]),
        SphericalPoint(b[1], b[2])
    )
end

arclength(s::SphericalSegment)::Float64 = arclength(s.a, s.b)

sphdist(s::SphericalSegment, R::Real=♂ᵣ)::Float64 = R*arclength(s)

function unitnormal(s::SphericalSegment)::SVector{3,Float64}
    unitnormal(sph2cart(s.a), sph2cart(s.b))
end

function checksegment(s::SphericalSegment)::Nothing
    checkpoint(s.a)
    checkpoint(s.b)
end

#==============================================================================
This stuff is for rotating points/segments over the sphere
==============================================================================#

export setuprotation, rotate, unrotate

function setuprotation(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T, tol::Real=1e-10) where {T<:Real}
    #coordinate integrity
    checkcoord(θ₁, ϕ₁)
    checkcoord(θ₂, ϕ₂)
    #account for cases with zero rotation to prevent NaN
    check₁ = isapprox(θ₁, θ₂, rtol=tol) & isapprox(ϕ₁, ϕ₂, rtol=tol)
    check₂ = (θ₁ < tol) & (θ₂ < tol)
    check₃ = (π - θ₁ < tol) & (π - θ₂ < tol)
    if check₁ | check₂ | check₃
        #zero rotation vector
        return (
            one(T),
            zero(T),
            unit(SVector{3,T}(one(T), one(T), one(T)))
        )
    end
    #cartesian vectors
    c₁ = sph2cart(θ₁, ϕ₁)
    c₂ = sph2cart(θ₂, ϕ₂)
    #angle between pole coordinates
    d = c₁ ⋅ c₂
    ψ = acos(c₁ ⋅ c₂)
    #axis of rotation
    k = unit(c₁ × c₂)
    return d, ψ, k
end

function rotate(θ::Float64,
                ϕ::Float64,
                d::T,
                ψ::T,
                k::SVector{3,T}
                )::NTuple{2,Float64} where {T<:Real}
    #cartesian location of rotating point
    v = sph2cart(θ, ϕ)
    #rotate
    s, c = sincos(ψ)
    w = v*c + (k × v)*s + k*(k ⋅ v)*(1.0 - d)
    #convert back to spherical coordinates, dropping radius
    cart2usph(w...)
end

function rotate(p::SphericalPoint, d, ψ, k)::SphericalPoint
    SphericalPoint(rotate(p.θ, p.ϕ, d, ψ, k))
end

function rotate(s::SphericalSegment, d, ψ, k)::SphericalSegment
    #just rotate both points
    SphericalSegment(
        rotate(s.a.θ, s.a.ϕ, d, ψ, k),
        rotate(s.b.θ, s.b.ϕ, d, ψ, k)
    )
end

rotate(θ, ϕ, args::T) where {T<:Tuple} = rotate(θ, ϕ, args...)

rotate(x, args::T) where {T<:Tuple} = rotate(x, args...)

function unrotate(θ::Float64,
                  ϕ::Float64,
                  d::T,
                  ψ::T,
                  k::SVector{3,T}
                  )::NTuple{2,Float64} where {T<:Real}
    #to reverse the rotation, k changes sign
    rotate(θ, ϕ, d, ψ, -k)
end

function unrotate(p::SphericalPoint, d, ψ, k)::SphericalPoint
    SphericalPoint(unrotate(p.θ, p.ϕ, d, ψ, k))
end

function unrotate(s::SphericalSegment, d, ψ, k)::SphericalSegment
    SphericalSegment(
        unrotate(s.a, d, ψ, k),
        unrotate(s.b, d, ψ, k)
    )
end

unrotate(θ, ϕ, args::T) where {T<:Tuple} = unrotate(θ, ϕ, args...)

unrotate(x, args::T) where {T<:Tuple} = unrotate(x, args...)

#==============================================================================
This type sets up a parameterized equation for a great circle through two
points, with periodic parameter range ∈ [0,2π]
# https://math.stackexchange.com/questions/1783746/equation-of-a-great-circle-passing-through-two-points
==============================================================================#

export GreatCircle
export sph, cart, colat
export 𝒪, 𝒪′, 𝒪′′, argmincolat
export 𝒟, 𝒟′, 𝒟′′, argmindist

struct GreatCircle
    v::SVector{3,Float64}
    w::SVector{3,Float64}
end

function GreatCircle(θ₁::Float64,
                     ϕ₁::Float64,
                     θ₂::Float64,
                     ϕ₂::Float64)
    @assert !isapprox(θ₁, θ₂, rtol=1e-13)
    @assert !isapprox(ϕ₁, ϕ₂, rtol=1e-13)
    v₁ = sph2cart(θ₁, ϕ₁)
    v₂ = sph2cart(θ₂, ϕ₂)
    d = v₁ ⋅ v₂
    f = sqrt(1 - d^2)
    α = -d/f
    β = 1/f
    w = α*v₁ + β*v₂
    GreatCircle(v₁, w)
end

function GreatCircle(a::SphericalPoint, b::SphericalPoint)
    GreatCircle(a.θ, a.ϕ, b.θ, b.ϕ)
end

GreatCircle(s::SphericalSegment) = GreatCircle(s.a, s.b)

#t is a parameter ∈ [0,2π] defining the circle
function (C::GreatCircle)(t)
    s, c = sincos(t)
    c*C.v + s*C.w
end

function sph(C::GreatCircle, t)::SphericalPoint
    SphericalPoint(cart2usph(C(t)...))
end

function xy(C::GreatCircle, t)
    x, y, _ = C(t)
    return x, y
end

function colat(C::GreatCircle, t)
    x, y, z = C(t)
    θ, _ = cart2sph(x, y, z)
    return θ
end

#objective function for size of vector projected on x-y plane
𝒪(C::GreatCircle, t) = sum(xy(C, t).^2)

𝒪′(C::GreatCircle, t) = derivative(x->𝒪(C,x), t)

𝒪′′(C::GreatCircle, t) = derivative(x->𝒪′(C,x), t)

#===
Find the point that minimizes x-y length using Newton's method.
Seems like this could be done analytically because the length
of the x-y projected great circle vector is a sine+cosine. This
is fast enough though, and quite accurate.
===#
function argmincolat(C::GreatCircle,
                     tol::Float64=1e-8,
                     maxiter::Int=100)::Float64
    #we know that the curve is sinusoidal   
    t = abs(𝒪′(C, 0.0)) < abs(𝒪′(C, π/4)) ? 0.0 : π/4
    Δ = Inf
    f′ = 𝒪′(C, t)
    f′′ = 𝒪′′(C, t)
    n::Int64 = 0
    while (abs(Δ) > tol) | (abs(f′) > tol)
        Δ = f′/f′′
        t = (t - Δ + 𝛕) % 𝛕
        f′ = 𝒪′(C, t)
        f′′ = 𝒪′′(C, t)
        n += 1
        n == maxiter && error("$maxiter iterations encountered $C")
    end
    #minimizer may be any of four correct solution
    T = (t, (t + π/2) % 𝛕, (t + π) % 𝛕, (t + 3π/2) % 𝛕)
    θ = map(t->colat(C,t), T)
    return T[argmin(θ)]
end

#objective function for distance between point and great circle
𝒟(C::GreatCircle, t, v::SVector{3,Float64}) = sum((C(t) - v).^2)

𝒟′(C::GreatCircle, t, v::SVector{3,Float64}) = derivative(x->𝒟(C,x,v), t)

𝒟′′(C::GreatCircle, t, v::SVector{3,Float64}) = derivative(x->𝒟′(C,x,v), t)

function argmindist(C::GreatCircle,
                    p::SphericalPoint,
                    tol::Float64=1e-8,
                    maxiter::Int=1000)::Float64
    #cartesian location of target point
    v = sph2cart(p)
    #newton's minimization
    t = 0.0
    Δ = Inf
    f′ = 𝒟′(C, t, v)
    f′′ = 𝒟′′(C, t, v)
    n::Int64 = 0
    while (abs(Δ) > tol) | (abs(f′) > tol)
        Δ = f′/f′′
        t = (t - Δ + 𝛕) % 𝛕
        f′ = 𝒟′(C, t, v)
        f′′ = 𝒟′′(C, t, v)
        n += 1
        n == maxiter && error("$maxiter iterations encountered $C")
    end
    𝒟(C, t, v) < 𝒟(C, t+π, v) ? t : (t + π) % 𝛕
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
    impacts::Int64 #number of registered impacts
    survived::Float64 #fraction of shoreline that survived
    destroyed::Float64 #fraction of shoreline destroyed
    segments::Vector{T} #surviving shoreline segments
    impactors::Set{Crater} #all craters registered as impacting the line
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
function segmentlengths(S::Vector{NTuple{2,Float64}},
                        θₛ::Float64, #segment latitude
                        R::Float64=♂ᵣ #sphere radius
                        )::Vector{Float64}

    #returned segment lengths
    seglen = Float64[]

    if length(S) == 1
        #a single segment should be a complete circle
        @assert S[1] == (0.0,𝛕)
        push!(seglen, 𝛕)
    else
        #multiple segments present
        for i ∈ 1:length(S)-1
            #segment length in radians
            push!(seglen, S[i][2] - S[i][1])
        end
        #final segment in radians
        Δϕ = S[end][2] - S[end][1]
        #check if it is distinct or wraps into the first seg
        if (S[1][1] == 0) & (S[end][2] == 𝛕)
            seglen[1] += Δϕ
        else
            push!(seglen, Δϕ)
        end
    end
    #convert to meters
    seglen .*= R*sin(θₛ)

    return seglen
end

#--------------------------------------
#iso-latitude representative shoreline

function root(crater::Crater,
              θ::Float64,
              Δϕ₁::Float64,
              Δϕ₂::Float64,
              tol::Float64=1e-8,
              maxiter::Int64=1000)::Float64
    d₁ = sphdist(c, θ, c.ϕ + Δϕ₁) - c.r
    d₂ = sphdist(c, θ, c.ϕ + Δϕ₂) - c.r
    Δϕ = Inf
    δϕ = Inf
    d = Inf
    n::Int64 = 0
    #secant method with stringent termination tolerance
    while (abs(d₂ - d₁) > tol) & (abs(δϕ) > tol)
        #approximate root
        δϕ = d₁*(Δϕ₂ - Δϕ₁)/(d₂ - d₁)
        Δϕ = Δϕ₁ - δϕ
        d = sphdist(c, θ, c.ϕ + Δϕ) - c.r
        #swaps
        Δϕ₁ = Δϕ₂
        Δϕ₂ = Δϕ
        d₁ = d₂
        d₂ = d
        #break on non-convergence
        n += 1
        n == maxiter && error("$maxiter iterations encountered, Δϕ₁=$Δϕ₁, Δϕ₂=$Δϕ₂, δϕ=$δϕ, d₁=$d₁, d₂=$d₂, d=$d, crater=$crater")
    end
    return Δϕ
end

function intersection(crater::Crater, θₛ::Real, R::Real=♂ᵣ)
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

function clip!(S::Vector{NTuple{2,T}}, sₙ::T, eₙ::T)::Nothing where {T<:Real}
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
                         rₑ::Float64, #ejecta scaling of radius
                         Δ::Float64 #required overlap distance for impact to register
                         )::SimulationResult
    #check coordinate boundaries
    @assert 0.0 <= θₛ <= π "shoreline colatitude (θₛ) must be ∈ [0,π]"
    #start a shoreline to take bites out of
    segments = NTuple{2,Float64}[(0.0,𝛕)]
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
            #register the impact
            push!(impactors, crater)
            #find the longitude intersection interval
            ϕ₁, ϕ₂ = intersection(crater, θₛ)
            if (ϕ₂ < ϕ₁)# & (abs(ϕ₁ - ϕ₂) > π/6)
                clip!(segments, 0., min(ϕ₁, ϕ₂))
                clip!(segments, max(ϕ₁, ϕ₂), 𝛕)
            else
                clip!(segments, ϕ₁, ϕ₂)
            end
        end
    end
    #compute the fraction surviving
    f = sum(seg -> seg[2] - seg[1], segments)/𝛕
    #construct the final result
    SimulationResult(
        length(impactors),
        f,
        1 - f,
        segments,
        impactors
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
                      latname::String="lat"
                      )::Vector{SphericalSegment}
    @assert minarc >= 0.0
    #read the table
    df = CSV.read(fn, DataFrame)
    #convert to radians
    θ, ϕ = latlon2sph(df[!,latname], df[!,lonname])
    N = length(θ)
    L = N - 1
    #accumulate segments
    S = SphericalSegment[]
    i = 1
    while i <= L
        #accumulate distance until exceeding mindist
        d = 0.0
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

function colatituderange(S::Vector{SphericalSegment})::NTuple{2,Float64}
    θa = map(s->s.a.θ, S)
    θb = map(s->s.b.θ, S)
    θmin = min(minimum(θa), minimum(θb))
    θmax = max(maximum(θa), maximum(θb))
    return θmin, θmax
end

function intersection(C::GreatCircle,
                      θ::Float64, #target colatitude
                      tₘ::Float64, #parameter of minimum colatitude
                      t::Float64, #other starting point
                      maxiter::Int=1000,
                      tol::Float64=1e-8)::Float64
    x₁ = tₘ
    x₂ = t
    y₁ = colat(C, x₁) - θ
    y₂ = colat(C, x₂) - θ
    @assert y₁*y₂ <= 0.0 "GreatCircle intersection not bounded by input parameters"
    @multiassign x, Δ, y = Inf
    n::Int64 = 0
    #regula falsi guarantees the solution is in the original interval
    while (abs(y) > tol) | (abs(Δ) > tol) 
        #approximate root
        Δ = y₁*(x₂ - x₁)/(y₂ - y₁)
        x = x₁ - Δ
        y = colat(C, x) - θ
        #reduce interval
        if y₁*y > 0
            x₁, y₁ = x, y
        else
            x₂, y₂ = x, y
        end
        #count
        n += 1
        if n == maxiter
            println("tₘ = $tₘ")
            println("θ = $θ")
            println("C = $C")
            error("$maxiter iterations encountered")
        end
    end
    return x
end

#parameter arguments where C intersects colatitude aᵣ
function intersections(C::GreatCircle, θ::Float64)::NTuple{2,Float64}
    #minimize the colatitude along the great circle
    tₘ = argmincolat(C)
    #first point where great circle meets θ
    t₁ = intersection(C, θ, tₘ, tₘ - π/6)
    #the other point where they meet, by symmetry
    t₂ = tₘ + (tₘ - t₁)
    #t₂ may be > 2π, but this is handled later
    return t₁, t₂
end

#function order(s::SphericalSegment, p::SphericalPoint)::Int64
#    𝓁a = arclength(s.a, p)
#    𝓁b = arclength(s.b, p)
#    println(𝓁a, ' ', 𝓁b, ' ', arclength(s) - (𝓁a + 𝓁b))
#    if isapprox(arclength(s), 𝓁a + 𝓁b, rtol=1e-10)
#        #point is inside the interval
#        return 0
#    else
#        if 𝓁a < 𝓁b
#            #point is on the a side
#            return -1
#        else
#            #point is on the b side
#            return 1
#        end
#    end
#    println(s)
#    println(p)
#    error("failure!")
#end
#
#function iscloser(x::SphericalPoint, a::SphericalPoint, b::SphericalPoint)::Bool
#    arclength(x, a) < arclength(x, b) ? true : false
#end

#function checkseg!(S::Vector{SphericalSegment},
#                   L::Int64,
#                   i::Int64)::Int64
#    if arclength(S[i]) < 1e-8
#        deleteat!(S, i)
#        return L - 1
#    end
#    return L
#end
#
#function updateseg!(S::Vector{SphericalSegment},
#                    L::Int64,
#                    i::Int64,
#                    a::SphericalPoint,
#                    b::SphericalPoint)::Int64
#    s = SphericalSegment(a, b)
#    if arclength(s) > 1e-8
#        S[i] = s
#        return L
#    else
#        deleteat!(S, i)
#        return L - 1
#    end
#end
#
#function clip!(S::Vector{SphericalSegment},
#               L::Int64, #original length of S
#               i::Int64, #index of segment under consideration
#               a::SphericalPoint,
#               b::SphericalPoint)::Int64
#    #====
#    All the points lie on the same great circle (or very close).
#    Strategy is to order the points using arclengths then check
#    for overlaps.
#    ====#
#    s = S[i]
#    o = order(s, a)
#    p = order(s, b)
#    #various overlap cases
#    if (o == 0) & (p == 0) #overlap inside the segment, break it
#        #see which point is closer to segment points
#        if iscloser(s.a, a, b)
#            println(1)
#            S[i] = SphericalSegment(s.a, a)
#            insert!(S, i+1, SphericalSegment(b, s.b))
#        else
#            println(2)
#            S[i] = SphericalSegment(s.a, b)
#            insert!(S, i+1, SphericalSegment(a, s.b))
#        end
#        L += 1
#        L = checkseg!(S, L, i)
#        L = checkseg!(S, L, i+1)
#    elseif ((o == -1) & (p == 1)) | ((o == 1) & (p == -1))
#        println(3)
#        #total destruction
#        deleteat!(S, i)
#        L -= 1
#    elseif (o == -1) & (p == 0)
#        println(4)
#        L = updateseg!(S, L, i, b, s.b)
#    elseif (o == 0) & (p == 1)
#        println(5)
#        L = updateseg!(S, L, i, s.a, a)
#    elseif (o == 0) & (p == -1)
#        println(6)
#        L = updateseg!(S, L, i, a, s.b)
#    elseif (o == 1) & (p == 0)
#        println(7)
#        L = updateseg!(S, L, i, s.a, b)
#    end
#    if arclength(S[i]) < 1e-12
#        println("$o $p")
#        println("$s\n$a\n$b")
#    end
#    return L
#end

function clip!(S::Vector{SphericalSegment},
               i::Int64,
               L::Int64,
               C::GreatCircle,
               sₙ::Float64,
               eₙ::Float64)#::Int64
    #segment under consideration
    s = S[i]
    #we know the t argument of the segment's points
    s, e = 0.0, arclength(s)
    #check various overlap cases
    if (s < sₙ) & (eₙ < e)
        
    elseif (sₙ <= s) & (e <= eₙ)
        
    elseif (sₙ <= s) & (s < eₙ < e)

    elseif (s < sₙ < e) & (e <= eₙ)

    end

end

function simulateimpacts(population::GlobalPopulation,
                         segments::Vector{SphericalSegment},
                         rₑ::Float64=1.0,
                         Δ::Float64=0.0
                         )::SimulationResult
    #check over segment coordinates
    for s ∈ segments
        checksegment(s)
    end
    L = length(segments)
    #find latitude range of segments
    θrange = colatituderange(segments)
    #make a copy of the segments before taking bites out of them
    segments = deepcopy(segments)
    #store craters that impact
    impactors = Set{Crater}()
    #iterate through the entire crater population
    for (count, crater) ∈ enumerate(population)
        #adjust crater radius for ejecta
        crater *= rₑ
        #short parameter names
        @unpack θ, ϕ, r = crater
        #arclength of crater radius
        aᵣ = r/♂ᵣ
        #============================================================
        The first check for intersection is simply whether the crater
        is so far from the colatitude range of the putative shoreline
        segments that it's impossible for it to touch any of them
        ============================================================#
        #if (θrange[1] - aᵣ) <= θ <= (θrange[2] + aᵣ)
            #set up rotation of this crater to the north pole
            rotation = setuprotation(θ, ϕ, 0.0, 0.0)
            #every segment has to be checked sadly
            i = 1
            while i <= L
                #========================================================
                An easy second check is whether the distance/arclength
                between crater center and segment endpoints far exceeds
                the crater's radius
                THIS HAS PROBLEMS???
                ========================================================#
                #𝓁a = arclength(s.a.θ, s.a.ϕ, θ, ϕ)
                #𝓁b = arclength(s.b.θ, s.b.ϕ, θ, ϕ)
                #if (𝓁a - aᵣ > π/4) & (𝓁b - aᵣ > π/4)
                    #rotate to put crater center at the north pole
                    x = rotate(segments[i], rotation)
                    #unit vector normal to plane through segment and origin
                    n = unitnormal(x)
                    #closest angle between plane and z axis
                    α = asin(abs(n[3]))
                    #the arclength of the overlap buffer
                    Δᵣ = Δ/♂ᵣ
                    #====================================================
                    The third check is whether the great circle formed
                    by the rotated segment runs through the crater,
                    which is now at the north pole. This is relatively
                    easy to check and should save a fair amount of time.
                    ====================================================#
                    if α + Δᵣ <= aᵣ
                        #================================================
                        By this stage optimization doesn't matter much 
                        because the bulk of the work is done rejecting
                        intersections before this branch is reached.
                        Things still need to be robust, of course.
                        ================================================#
                        #paramaterize the rotated segment's great circle
                        C = GreatCircle(x)
                        #parameter values where C intersects the crater
                        sₙ, eₙ = intersections(C, aᵣ)
                        #we know the param values of segment points immediately
                        s, e = 0.0, arclength(x)
                        println("$s $e $sₙ $eₙ")
                        

                        push!(impactors, crater)
                    end
                #end
                i += 1
            end
        #end
        #occasionally refresh the colatitude range
        if count % 10_000 == 0
            θrange = colatituderange(segments)
        end
    end
    SimulationResult(length(impactors), NaN, NaN, segments, impactors)
end

function simulateimpacts(population::GlobalPopulation,
                         segments::Vector{SphericalSegment},
                         rₑ::Real,
                         Δ::Real)::SimulationResult
    simulateimpacts(population, segments, Float64(rₑ), Float64(Δ))
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
