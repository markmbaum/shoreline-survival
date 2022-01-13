module ShorelineSurvival

import Base.*, Base.big
using Base.Threads: @threads, nthreads
using LinearAlgebra: ⋅, ×
using PrettyTables
using UnPack
using Random: AbstractRNG, MersenneTwister, seed!
using StaticArrays
using MultiAssign
using Formatting
using CSV
using DataFrames

include("geometry.jl")
include("craters.jl")
include("utilities.jl")

#==============================================================================
The following functions and definitions handle the simulation of craters
impacting a hypothetical shoreline on the global scale
==============================================================================#

export GlobalResult
export survived, destroyed, segdistances, gapdistances

struct GlobalResult{T}
    A₀::Float64 #original total arclength of segments
    A::Float64 #final total arclength of segments
    impacts::Int64 #number of registered impacts
    segments::Vector{T} #surviving shoreline segments
    impactors::Vector{Crater} #all craters registered as impacting the line
end

function Base.show(io::IO, res::GlobalResult{T}) where {T}
    println(io, "GlobalResult{$T}")
    println(io, "  $(res.impacts) impacts registered")
    A₀ = round(res.A₀, sigdigits=6)
    println(io, "  initial Σarclength = $A₀ radians")
    A = round(res.A, sigdigits=6)
    println(io, "  final   ∑arclength = $A radians")
    f = round(100*survived(res), sigdigits=8)
    println(io, "  $f % survived")
    f = round(100*destroyed(res), sigdigits=8)
    print(io, "  $f % destroyed")
end

survived(res::GlobalResult) = res.A/res.A₀

destroyed(res::GlobalResult) = 1 - survived(res)

#computes segment lengths (in meters) of an impacted shoreline
function segdistances(S::Vector{NTuple{2,Float64}},
                      θₛ::Float64, #segment latitude
                      R::Float64=♂ᵣ #sphere radius
                      )::Vector{Float64}
    if length(S) == 0
        a = [0.0]
    elseif length(S) == 1
        a = [S[1][2] - S[1][1]]
    else
        a = map(s->s[2]-s[1], S)
        #check if first and last segments actually wrap
        if (S[1][1] == 0) & (S[end][2] == τ)
            a[1] += pop!(a)
        end
    end
    #scale by radius and latitude to get distance in meters
    return R*sin(θₛ)*a
end

function segdistances(segments::Vector{SphericalSegment{T}}, R::Float64=♂ᵣ) where {T}
    #escape hatch the empty case
    length(segments) == 0 && return [zero(T)]
    #use big arithmetic for extra arclength accuracy
    S = big.(segments)
    #assume the segments are in order
    a = ∠.(S)
    𝓁 = T[a[1]] #start a list of segment arclengths
    for i ∈ 2:length(S)
        if commonendpoint(S[i], S[i-1])
            𝓁[end] += a[i]
        else
            push!(𝓁, convert(T,a[i]))
        end
    end
    #handle possible wrapping
    if commonendpoint(S[1], S[end]) & (length(𝓁) > 1)
        𝓁[1] += pop!(𝓁)
    end
    #remember to apply the radius
    return R*𝓁
end

function segdistances(res::GlobalResult{T}, args...) where {T}
    segdistances(res.segments, args...)
end

#computes gap lengths (in meters) of an impacted shoreline
function gapdistances(S::Vector{NTuple{2,Float64}},
                      θₛ::Float64, #segment latitude
                      R::Float64=♂ᵣ #sphere radius
                      )::Vector{Float64}
    if length(S) == 0
        g = [τ]
    elseif length(S) == 1
        g = [τ - (S[1][2] - S[1][1])]
    else
        g = Float64[]
        for i ∈ 2:length(S)
            e = S[i-1][2]
            s = S[i][1]
            if e != s
                #segments are not connected, store the gap
                push!(g, s - e)
            end
        end
        #check for lack of wrapping
        s = S[1][1]
        e = S[end][2]
        if (s != 0) & (e != τ)
            push!(g, s + τ - e)
        end
    end
    #scale by radius and latitude to get distance in meters
    return g*R*sin(θₛ)
end

#====
This function may be unrepresentative for very strange cases, like a single
small segment floating alone on the sphere. It assumes minimum distance
between non-overlapping segments, so cases where the gap arc length should be
greater than π will be miscalculated. These cases should be very unlikely,
however.
====#
function gapdistances(segments::Vector{SphericalSegment{T}}, R::Float64=♂ᵣ) where {T}
    #escape hatch the empty case
    length(segments) == 0 && [convert(T,NaN)]
    #use big arithmetic
    S = big.(segments)
    #assume the segments are in order
    g = T[]
    for i ∈ 2:length(S)
        if !commonendpoint(S[i], S[i-1])
            #segments are not connected, store the gap
            push!(g, ∠(S[i-1].b, S[i].a))
        end
    end
    #check for lack of wrapping
    if !commonendpoint(S[end], S[1])
        push!(g, ∠(S[end].b, S[1].a))
    end
    #remember to apply the radius
    return R*g
end

function gapdistances(res::GlobalResult{T}, args...) where {T}
    gapdistances(res.segments, args...)
end

#--------------------------------------
#iso-latitude representative shoreline

export intersection

function intersection(θ::T, ϕ::T, r::T, θₛ::T, R::T)::Float64 where {T<:Real}
    #equal to dot product of crater center vector and solution pt vector
    C = cos(r/R)
    #magnitude of any x-y vector in the θₛ circle
    S = sin(θₛ)
    #height of the θₛ circle
    z = cos(θₛ)
    #cartesian circle center
    a, b, c = sph2cart(θ, ϕ)
    #precompute some things
    a², b², c² = a^2, b^2, c^2
    C², S², z² = C^2, S^2, z^2
    #x coordinate of one point (other point flips sign of sqrt term)
    x = (a*C - a*c*z - sqrt(-b²*C² + a²*b²*S² + b^4*S² + 2*b²*c*C*z - b²*c²*z²))/(a² + b²)
    #y coordinates are the same for both points
    y = (C - (a²*C)/(a² + b²) - c*z + (a²*c*z)/(a² + b²) + (a*sqrt(-b²*(C² - a²*S² - b²*S² - 2*c*C*z + c²*z²)))/(a² + b²))/b
    #convert back to a longitude interval
    @inbounds Float64(abs(cart2sph(x, y, z)[2] - ϕ))
end

function intersection(crater::Crater{Float64}, θₛ::Float64, R::Float64)::NTuple{2,Float64}
    @unpack θ, ϕ, r = crater
    #slower to use extended precision but very little impact on overall speed
    Δϕ = intersection(θ, ϕ, r, θₛ, R)
    #create a longitude interval with values ∈ [0,2π]
    return ↻(ϕ - Δϕ), ↻(ϕ + Δϕ)
end

function clip!(S::Vector{NTuple{2,T}}, sₙ::T, eₙ::T)::Bool where {T<:Real}
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

function globalsimulation(population::GlobalPopulation,
                          θₛ::Float64, #shoreline colatitude [0,π]
                          rₑ::Float64, #ejecta scaling of radius
                          Δ::Float64 #required overlap distance for impact to register
                          )::GlobalResult
    #check coordinate boundaries
    @assert 0.0 <= θₛ <= π "shoreline colatitude (θₛ) must be ∈ [0,π]"
    #start a shoreline to take bites out of
    segs = NTuple{2,Float64}[(0.0,τ)]
    #store craters that impact
    impactors = Set{Crater}()
    #now go through each crater, chopping up the shoreline as necessary
    for crater ∈ population
        #adjust radius for ejecta and unpack (must be done by creating new Crater!)
        crater *= rₑ
        @unpack θ, ϕ, r = crater
        #distance from crater center to line
        dₛ = ♂ᵣ*abs(θₛ - θ)
        #check if the crater overlaps the line enough
        if dₛ < r - Δ
            #find the longitude intersection interval
            ϕ₁, ϕ₂ = intersection(crater, θₛ, ♂ᵣ)
            #clip overlapping portions
            if ϕ₂ < ϕ₁ #intersection interval wraps over 2π
                if clip!(segs, 0., min(ϕ₁, ϕ₂)) | clip!(segs, max(ϕ₁, ϕ₂), τ)
                    push!(impactors, crater)
                end
            else
                clip!(segs, ϕ₁, ϕ₂) && push!(impactors, crater)
            end
        end
    end
    #initial arclenth of isolatitude ring
    A₀ = τ*sin(θₛ)
    #total arclength of surviving segments
    A = sin(θₛ)*(length(segs) >= 1 ? sum(x->x[2]-x[1], segs) : 0.0)
    #construct the final result
    GlobalResult(A₀, A, length(impactors), segs, collect(impactors))
end

#barrier function
function globalsimulation(population::GlobalPopulation, θₛ::Real, rₑ::Real, Δ::Real)
    globalsimulation(population, Float64(θₛ), Float64(rₑ), Float64(Δ))
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
        #accumulate distance until exceeding minimum arc length
        d = zero(T)
        j = i
        while (d <= minarc) & (j < N)
            d = ∠(θ[i], ϕ[i], θ[j+1], ϕ[j+1])
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

function newseg(𝓋₁::SVector{3,T},
                𝓋′::SVector{3,T},
                t₁::Float64,
                t₂::Float64)::CartesianSegment where {T<:Real}
    CartesianSegment(
        𝓋₁*cos(t₁) + 𝓋′*sin(t₁), #evaluates great circle at t₁
        𝓋₁*cos(t₂) + 𝓋′*sin(t₂)
    )
end

function clip!(csegs::Vector{CartesianSegment{Float64}},
               𝓊::Vector{SVector{3,Float64}},
               i::Int64,
               𝓋₁::SVector{3,Float64},
               𝓋′::SVector{3,Float64},
               𝓁ᵢ::Float64,
               sₙ::Float64,
               eₙ::Float64)::Tuple{Int64,Bool}
    #always have the parameter arguments of the segment's end points
    s, e = 0.0, 𝓁ᵢ
    #check various overlap cases
    ΔL = 0
    case = overlapcase(s, e, sₙ, eₙ)
    if case == 1
        #split the segment, first part
        csegs[i] = newseg(𝓋₁, 𝓋′, s, sₙ)
        #second part
        insert!(csegs, i+1, newseg(𝓋₁, 𝓋′, eₙ, e))
        insert!(𝓊, i+1, 𝓊[i])
        ΔL = 1
    elseif case == 2
        #intersection contains the segment, discard the seg
        deleteat!(csegs, i)
        deleteat!(𝓊, i)
        ΔL = -1
    elseif case == 3
        #overlap on the lower side, crop up
        csegs[i] = newseg(𝓋₁, 𝓋′, eₙ, e)
    elseif case == 4
        #overlap on the upper side, crop down
        csegs[i] = newseg(𝓋₁, 𝓋′, s, sₙ)
    end
    return ΔL, (case == 0) ? false : true
end

function subglobalsimulation(population::GlobalPopulation,
                             segs::Vector{SphericalSegment{𝒯}},
                             rₑ::Float64=1.0,
                             Δ::Float64=0.0,
                             minarc::Float64=1/♂ᵣ) where {𝒯<:Real}
    #check over segment coordinates
    foreach(checksegment, segs)
    #initial number of segments
    L = length(segs)
    #find latitude range of segments
    θmin, θmax = colatrange(segs)
    #store initial sum of segment arclengths to compare with at the end
    A₀ = sum(map(∠, segs))
    #keep a cartesian mirror of the segments to speed up first filter
    csegs::Vector{CartesianSegment{𝒯}} = map(CartesianSegment, segs)
    #pre-compute unit vectors normal to the original segments
    𝓊::Vector{SVector{3,𝒯}} = map(unitnormal, csegs)
    #store craters that impact
    impactors = Set{Crater}()
    #the arclength of the overlap buffer
    Δᵣ = Δ/♂ᵣ
    #iterate through the entire crater population
    for crater ∈ population
        #adjust radius for ejecta and unpack (must be done by creating new Crater!)
        crater *= rₑ
        @unpack θ, ϕ, r = crater
        #cartesian crater center
        χ = sph2cart(θ, ϕ)
        #arclength of crater radius
        𝓁ᵣ = r/♂ᵣ
        #store for check on plane intersection, preventing lots of asin evals
        𝒮⁺ = sin(𝓁ᵣ - Δᵣ)
        𝒮⁻ = -𝒮⁺
        #============================================================
        The first check for intersection is simply whether the crater
        is so far from the colatitude range of the line segments that
        it's impossible for it to touch any of them
        ============================================================#
        if (θmin - 𝓁ᵣ) <= θ <= (θmax + 𝓁ᵣ)
            #sadly, every segment now has to be checked
            i = 1
            while i <= L
                #====================================================
                The second check is whether the plane formed
                by the rotated segment runs through the crater,
                which is now at the north pole. This is relatively
                easy to check and should save a fair amount of time.
                This is simultaneously a check that the overlap meets
                the minimum requirement Δ.
                ====================================================#
                if @inbounds 𝒮⁻ < 𝓊[i] ⋅ χ < 𝒮⁺
                    #========================================================
                    The third check is whether the distance/arclength
                    between crater center and segment endpoints far exceeds
                    the crater's radius
                    ========================================================#
                    @inbounds cᵢ = csegs[i]
                    if (∠(χ, cᵢ.a) + 𝓁ᵣ < π/4) & (∠(χ, cᵢ.b) + 𝓁ᵣ < π/4)
                        #================================================
                        By this stage optimization doesn't matter much 
                        because the bulk of the work is done rejecting
                        intersections before this branch is reached.
                        Things still need to be robust, of course.
                        ================================================#
                        #arclength of the actual segment
                        𝓁ᵢ = @inbounds ∠(cᵢ)
                        #check if the segment is too small to keep
                        if 𝓁ᵢ < minarc
                            deleteat!(csegs, i)
                            deleteat!(𝓊, i)
                            L -= 1
                            i -= 1
                        else
                            #parameter values where C intersects the crater
                            # https://math.stackexchange.com/questions/4330547/intersection-of-circle-and-geodesic-segment-on-sphere
                            𝓋₁ = cᵢ.a
                            𝓋₂ = cᵢ.b
                            𝓋′ = unit(𝓋₂ - 𝓋₁*(𝓋₁ ⋅ 𝓋₂))
                            A = χ ⋅ 𝓋₁
                            B = χ ⋅ 𝓋′
                            t₀ = atan(B, A)
                            α = cos(𝓁ᵣ)/sqrt(A^2 + B^2)
                            Δt = acos(α > 1.0 ? 1.0 : α) #there might be value *barely* greater than 1
                            t₁ = t₀ - Δt
                            t₂ = t₀ + Δt
                            #sanity check that intersection segment is not larger than crater
                            Δt > 2𝓁ᵣ && @assert (Δt - 2𝓁ᵣ)/(2𝓁ᵣ) < 1e-6
                            #now check for genuine overlap
                            ΔL, impacted = clip!(csegs, 𝓊, i, 𝓋₁, 𝓋′, 𝓁ᵢ, t₁, t₂)
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
    end
    #convert final segments back to spherical coordinates
    segs::Vector{SphericalSegment{𝒯}} = map(SphericalSegment, csegs)
    #final sum of segment arclengths
    A = sum(map(∠, segs))
    #final construction
    GlobalResult(A₀, A, length(impactors), segs, collect(impactors))
end

function globalsimulation(population::GlobalPopulation,
                          segs::Vector{SphericalSegment{Float64}},
                          rₑ::Float64=1.0,
                          Δ::Float64=0.0,
                          minarc::Float64=1/♂ᵣ)
    #number of groups/threads
    N = nthreads()
    #split segments up into groups
    subsegs = makechunks(segs, N)
    #work on each group in parallel
    res = Vector{GlobalResult{SphericalSegment{Float64}}}(undef, N)
    #@threads 
    for i ∈ 1:N
        res[i] = subglobalsimulation(
            deepcopy(population),
            subsegs[i],
            rₑ,
            Δ,
            minarc
        )
    end
    #put all the subresults together
    GlobalResult(
        sum(getfield.(res, :A₀)),
        sum(getfield.(res, :A)),
        sum(getfield.(res, :impacts)),
        vcat(getfield.(res, :segments)...),
        vcat(getfield.(res, :impactors)...)
    )
end

#barrier function
function globalsimulation(population::GlobalPopulation,
                         segments::Vector{SphericalSegment{Float64}},
                         rₑ::Real,
                         Δ::Real)
    globalsimulation(population, segments, Float64(rₑ), Float64(Δ))
end

#--------------------------------------
#convenience function

export globalsimulation

function globalsimulation(t::Real, #time [Ga]
                          shoreline::𝒯, #putative shoreline segments or latitude
                          rₑ::Real=1.0, #ejecta scaling of radius
                          Δ::Real=0.0; #required overlap distance for impact to register
                          rmin::Real=1e3, #smallest allowed crater radius [m]
                          nmax::Real=1_000_000, #maximum craters in bins, default small value
                          seed=1,
                          show::Bool=false) where {𝒯}
    #check overlap distance
    @assert Δ >= 0 "overlap distance (Δ) must be positive"
    #check ejecta radius multiple
    @assert rₑ >= 1 "ejecta radius multiple (rₑ) must be greater than or equal to 1"
    #start up the crater population (impossible to have impacts where r < Δ)
    population = GlobalPopulation(t, rmin=max(rmin,Δ), nmax=nmax, seed=seed)
    #print the crater population table if desired
    show && println(population)
    #send the craters!
    globalsimulation(population, shoreline, rₑ, Δ)
end

end
