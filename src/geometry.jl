#==============================================================================
This file contains functions for doing various things in spherical geometry.
I use a colatitude coordinate θ (theta) ∈ [0,π] and a longitude
coordinate ϕ (phi) ∈ [0,2π].
==============================================================================#

const τ = 2π
export τ

#--------------------------------------
export sphrand

function sphrand(rng::AbstractRNG)::NTuple{2,Float64}
    θ = acos(1 - 2*rand(rng))
    ϕ = τ*rand(rng)
    return θ, ϕ
end

function sphrand()::NTuple{2,Float64}
    θ = acos(1 - 2*rand())
    ϕ = τ*rand()
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

export ∠, sphdist

#this is the arclength, assuming vectors have length 1
function ∠(c₁::SVector{3,T}, c₂::SVector{3,T}) where {T}
    d = c₁ ⋅ c₂
    if d > 1
        return zero(T)
    elseif d < -1
        return convert(T,π)
    end
    acos(d)
end

function ∠(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where {T}
    ∠(sph2cart(θ₁, ϕ₁), sph2cart(θ₂, ϕ₂))
end

sphdist(θ₁, ϕ₁, θ₂, ϕ₂, R) = R*∠(θ₁, ϕ₁, θ₂, ϕ₂)

#--------------------------------------
#misc

export ↻, unit, sphcirc, unitnormal

#wraps an angle into [0,2π] and appears to be quicker than using remainder
function ↻(θ)
    while θ < 0; θ += τ; end
    while θ > τ; θ -= τ; end
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
    @inbounds for (i,ψ) ∈ enumerate(LinRange(0, τ, N))
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
    @assert 0.0 <= ϕ <= τ
    nothing
end

#==============================================================================
Here two simple types for spherical geometry are defined, along with
some basic operations on them as wrappers of the functions above.
==============================================================================#

#--------------------------------------
export SphericalPoint

struct SphericalPoint{T<:AbstractFloat}
    θ::T
    ϕ::T
    function SphericalPoint{T}(θ::T, ϕ::T) where {T<:AbstractFloat}
        @assert 0 <= θ <= π "colatitude must be ∈ [0,π]"
        @assert 0 <= ϕ <= τ "longitude must be ∈ [0,2π]"
        new{T}(θ, ϕ)
    end
end

SphericalPoint(θ::T, ϕ::T) where {T<:AbstractFloat} = SphericalPoint{T}(θ,ϕ)

Base.show(io::IO, p::SphericalPoint) = print(io, "(θ=$(p.θ), ϕ=$(p.ϕ))")

SphericalPoint(x::NTuple{2,T}) where {T} = @inbounds SphericalPoint{T}(x[1], x[2])

sph2cart(p::SphericalPoint) = sph2cart(p.θ, p.ϕ)

∠(a::SphericalPoint{T}, b::SphericalPoint{T}) where {T} = ∠(a.θ, a.ϕ, b.θ, b.ϕ)

checkpoint(p::SphericalPoint)::Nothing = checkcoord(p.θ, p.ϕ)

Base.big(p::SphericalPoint{T}) where {T} = SphericalPoint(big(p.θ), big(p.ϕ))

#--------------------------------------
export SphericalSegment
export commonendpoint

struct SphericalSegment{T<:AbstractFloat}
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

function Base.show(io::IO, s::SphericalSegment{T}) where {T}
    print(io, "SphericalSegment{$T}\n  a = $(s.a)\n  b = $(s.b)")
end

∠(s::SphericalSegment) = ∠(s.a, s.b)

sphdist(s::SphericalSegment, R::Real=♂ᵣ) = ∠(s)*R

sph2cart(s::SphericalSegment) = sph2cart(s.a), sph2cart(s.b)

function checksegment(s::SphericalSegment, maxarc=π/6)::Nothing
    checkpoint(s.a)
    checkpoint(s.b)
    𝓁 = ∠(s)
    if 𝓁 > maxarc
        p = round(100*𝓁/τ, sigdigits=4)
        error("unusually large segment with arclength=$𝓁 or ~$p % of 2π")
    end
    nothing
end

function commonendpoint(s₁::SphericalSegment{T}, s₂::SphericalSegment{T})::Bool where {T}
    c₁ = (sph2cart(s₁.a), sph2cart(s₁.b))
    c₂ = (sph2cart(s₂.a), sph2cart(s₂.b))
    for p₁ ∈ c₁, p₂ ∈ c₂
        #a tiny bit of wiggle room on the unit sphere
        all(p₁ ≈ p₂) && return true
    end
    return false
end

Base.big(s::SphericalSegment{T}) where {T} = SphericalSegment(big(s.a), big(s.b))

#--------------------------------------
export CartesianSegment

struct CartesianSegment{T}
    a::SVector{3,T}
    b::SVector{3,T}
end

function CartesianSegment(s::SphericalSegment{T}) where {T}
    CartesianSegment{T}(sph2cart(s.a), sph2cart(s.b))
end

function SphericalSegment(c::CartesianSegment{T}) where {T}
    SphericalSegment(cart2usph(c.a), cart2usph(c.b))
end

∠(c::CartesianSegment) = ∠(c.a, c.b)

unitnormal(c::CartesianSegment) = unitnormal(c.a, c.b)

function minplanecolat(c::CartesianSegment)
    n = unitnormal(c.a, c.b)
    z = @inbounds abs(n[3])
    return asin(z)
end

Base.big(c::CartesianSegment{T}) where {T} = CartesianSegment(big.(c.a), big.(c.b))

#==============================================================================
This type sets up a parameterized equation for a great circle through two
points, with periodic parameter range ∈ [0,2π]
# https://math.stackexchange.com/questions/1783746/equation-of-a-great-circle-passing-through-two-points
==============================================================================#

export GreatCircle
export sph, colat

struct GreatCircle{T}
    v::SVector{3,T}
    w::SVector{3,T}
    function GreatCircle(v₁::SVector{3,T}, v₂::SVector{3,T}) where {T<:Real}
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

function GreatCircle(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where {T<:Real}
    GreatCircle(sph2cart(θ₁, ϕ₁), sph2cart(θ₂, ϕ₂))
end

function GreatCircle(a::SphericalPoint{T}, b::SphericalPoint{T}) where {T<:Real}
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
    θ, _ = cart2usph(x, y, z)
    return θ
end