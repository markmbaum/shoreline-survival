using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using DataFrames
using CSV
using PyPlot
using Statistics

pygui(true)

##

function fixedges!(S)::Nothing
    s = S[1]
    if s.a.ϕ > π
        S[1] = SphericalSegment(SphericalPoint(s.a.θ, 0.0), s.b)
    end
    s = S[end]
    if s.b.ϕ < π
        S[end] = SphericalSegment(s.a, SphericalPoint(s.b.θ, τ))
    end
    nothing
end

function meanline(df::DataFrame, t, rₑ)::Int64
    #filter
    sl = df[df.re .== rₑ,:]
    sl = sl[sl.t .== t,:]
    #mean survival frac
    f = mean(sl.survived)
    #index of result closest to the mean
    sl[argmin(@. abs(sl.survived - f)), :idx]
end

function rectangle!(ax, x₁, x₂)::Nothing
    y₁, y₂ = ax.get_ylim()
    ax.add_patch(
        matplotlib.patches.Rectangle(
            (x₁, y₁),
            x₂ - x₁,
            y₂ - y₁,
            zorder=-1,
            color="gray",
            ec=nothing,
            alpha=0.25
        )
    )
    nothing
end

function zoomlines!(ax₁, ax₂, x₁, x₂)::Nothing
    yl₁ = ax₁.get_ylim()
    ax₁.add_artist(
        matplotlib.patches.ConnectionPatch(
            (x₁, yl₁[1]),
            (0, 1),
            "data",
            "axes fraction",
            ax₁,
            ax₂,
            alpha=0.7
        )
    )
    ax₁.add_artist(
        matplotlib.patches.ConnectionPatch(
            (x₂, yl₁[1]),
            (1, 1),
            "data",
            "axes fraction",
            ax₁,
            ax₂,
            alpha=0.7
        )
    )
    nothing
end

ϕ2lon(ϕ) = ϕ*(180/π) - 180

θ2lat(θ) = -θ*(180/π) + 90

scale(ϕ₁, ϕ₂) = Int(round(♂ᵣ*(ϕ₂ - ϕ₁)/1e3))

## load the original shoreline coordiantes

S₀ = readsegments(
    datadir(
        "exp_pro",
        "parker_1989_contact_1a.csv"
    ),
    minarc=0.05
)
θ₀, ϕ₀ = flattensegments(S₀)

## load segments

dirseg = datadir("sims", "mapped-segments")
segments = Dict{Int64,Vector{SphericalSegment{Float64}}}()
paths = filter(p->!occursin(".txt", p), readdir(dirseg, join=true))
for path ∈ paths
    idx = parse(Int64, split(path, '_')[end])
    segments[idx] = loadsegments(path)
    fixedges!(segments[idx])
end

## results table

df = DataFrame(CSV.File(datadir("sims", "mapped.csv")))
df[!,:idx] = 1:size(df,1)
df = df[df.t .== 4,:]

##

fig, axs = subplots(4, 1, figsize=(4,5))
θ, ϕ = flattensegments(segments[meanline(df, 4, 1.5)])
color = plt.cm.cool(0.5)
z = [-1., 1.]
zoom₁ = π .+ (π/9)*z
zoom₂ = π .+ (π/90)*z
zoom₃ = π .+ (π/900)*z

axs[1].plot(ϕ₀, θ2lat.(θ₀), linewidth=0.5, color="k", alpha=0.7)
axs[1].plot(ϕ, θ2lat.(θ), linewidth=2.2, color=color)
axs[1].set_xlim(0, τ)
axs[1].set_xticks([0,τ])
axs[1].set_xticklabels(round.(ϕ2lon.((0,τ)), sigdigits=2), fontsize=7)
axs[1].set_title("global scale", fontsize=9)
rectangle!(axs[1], zoom₁[1], zoom₁[2])
zoomlines!(axs[1], axs[2], zoom₁...)

m = @. (zoom₁[1] <= ϕ <= zoom₁[2]) | isnan(ϕ)
axs[2].plot(ϕ[m], θ2lat.(θ[m]), linewidth=2.2, color=color)
axs[2].set_xlim(zoom₁)
axs[2].set_xticks([zoom₁[1], zoom₁[2]])
axs[2].set_xticklabels(round.(ϕ2lon.(zoom₁), sigdigits=2), fontsize=7)
axs[2].set_title("~$(scale(zoom₁...)) km", fontsize=9)
rectangle!(axs[2], zoom₂[1], zoom₂[2])
zoomlines!(axs[2], axs[3], zoom₂...)

m = @. (zoom₂[1] <= ϕ <= zoom₂[2]) | isnan(ϕ)
axs[3].plot(ϕ[m], θ2lat.(θ[m]), linewidth=2.2, color=color)
axs[3].set_xlim(zoom₂)
axs[3].set_xticks([zoom₂[1], zoom₂[2]])
axs[3].set_xticklabels(round.(ϕ2lon.(zoom₂), sigdigits=2), fontsize=7)
axs[3].set_title("~$(scale(zoom₂...)) km", fontsize=9)
rectangle!(axs[3], zoom₃[1], zoom₃[2])
zoomlines!(axs[3], axs[4], zoom₃...)

m = @. (zoom₃[1] <= ϕ <= zoom₃[2]) | isnan(ϕ)
axs[4].plot(ϕ[m], θ2lat.(θ[m]), linewidth=2.2, color=color)
axs[4].set_xlim(zoom₃)
axs[4].set_xticks([zoom₃[1], zoom₃[2]])
axs[4].set_xticklabels(round.(ϕ2lon.(zoom₃), sigdigits=2), fontsize=7)
axs[4].set_title("~$(scale(zoom₃...)) km", fontsize=9)
axs[4].set_xlabel("Longitude [deg]", fontsize=9)

for ax in axs
    ax.set_yticks([])
    ax.grid(false)
    for spine ∈ ("left", "right", "top", "bottom")
        ax.spines[spine].set_visible(true)
    end
end
#fig.supylabel("Latitude [deg]")
fig.tight_layout()
plt.subplots_adjust(hspace=0.45)
fig.savefig(plotsdir("results", "fractal_segments"), dpi=500)
plt.close(fig)