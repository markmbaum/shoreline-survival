#==============================================================================
This file contains various stand-alone functions with general applicability
==============================================================================#

function overlapcase(s::T, e::T, sₙ::T, eₙ::T)::Int64 where {T<:Real}
    if (sₙ >= e) | (eₙ <= s)
        #no overlap
        return 0
    elseif (s < sₙ) & (eₙ < e)
        #new interval is inside
        return 1
    elseif (sₙ <= s) & (e <= eₙ)
        #contained
        return 2
    elseif (sₙ <= s) & (s <= eₙ <= e)
        #overlap on the lower side
        return 3
    elseif (s <= sₙ <= e) & (e <= eₙ)
        #overlap on the upper side
        return 4
    else
        println("s=$s, e=$e, sₙ=$sₙ, eₙ=$eₙ")
        error("overlap case failure")
    end
end

function makechunks(X::AbstractVector{T}, n::Int) where {T}
    L = length(X)
    c = L ÷ n
    Y = Vector{Vector{T}}(undef, n)
    idx = 1
    for i ∈ 1:n-1
        Y[i] = X[idx:idx+c-1]
        idx += c 
    end
    Y[end] = X[idx:end]
    return Y
end