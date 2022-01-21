using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using Base.Threads: @threads
using Statistics: std

##

#area of a grid box rectangular in latitude and longitude
# colatitude θ ∈ [0,π]
# longitude ϕ ∈ [0,2π]
cellarea(Δϕ, θ₁, θ₂) = Δϕ*(cos(θ₁) - cos(θ₂))

function sphranddensity(npt::Int, nθ::Int, nϕ::Int)
    θ, ϕ = sphrand(npt)
    N = length(θ)
    S = zeros(Int64, nθ, nϕ)
    A = zeros(nθ, nϕ)
    Δθ = π/nθ
    Δϕ = 2π/nϕ
    @threads for i ∈ 1:nθ
        θ₁ = Δθ*(i-1)
        θ₂ = Δθ*i
        @inbounds for j ∈ 1:nϕ
            ϕ₁ = Δϕ*(j-1)
            ϕ₂ = Δϕ*j
            s::Int64 = 0
            for k ∈ 1:N
                θₖ = θ[k]
                ϕₖ = ϕ[k]
                if (θ₁ < θₖ < θ₂) && (ϕ₁ < ϕₖ < ϕ₂)
                    s += 1
                end
            end
            S[i,j] = s
            A[i,j] = cellarea(Δϕ, θ₁, θ₂)
        end
    end
    ρ = nθ*nϕ*S./(A*N)
    return ρ
end

##

npt = 1000
while npt < 1e9
    ρ = sphranddensity(npt, 10, 3)
    println("n = $npt\nstd = $(std(ρ))\n")
    npt *= 10
end