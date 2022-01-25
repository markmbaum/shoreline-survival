using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using DataFrames
using CSV
using PyPlot
using Statistics

pygui(true)

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
for path ∈ readdir(dirseg, join=true)
    idx = parse(Int64, split(path, '_')[end])
    segments[idx] = loadsegments(path)
end

## results table

df = DataFrame(CSV.File(datadir("sims", "mapped.csv")))
df[!,:idx] = 1:size(df,1)
df = df[df.t .== 4,:]

## plots!

fig = figure()
i = rand(1:size(df,1))
idx = df[i,:idx]
println(df[i,1:6])
θ, ϕ = flattensegments(segments[idx])
plot(ϕ, θ, color="C0", alpha=0.8, linewidth=1.5)

##

function meanline(df::DataFrame, t, rₑ)::Int64
    #filter
    sl = df[df.re .== rₑ,:]
    sl = sl[sl.t .== t,:]
    #mean survival frac
    f = mean(sl.survived)
    #index of result closest to the mean
    sl[argmin(@. abs(sl.survived - f)), :idx]
end

##

fig, axs = subplots(2, 1, figsize=(8,5))
urₑ = sort(unique(df.re))
colors = plt.cm.cool.(LinRange(0, 1, length(urₑ)))
zoomlim = (15π/16, 17π/16)
for (i,rₑ) ∈ enumerate(urₑ)
    θ, ϕ = flattensegments(segments[meanline(df, 4, rₑ)])
    axs[1][:plot](
        ϕ,
        -θ .+ 0.9*i,
        linewidth=1.5,
        color=colors[i],
        label=rₑ
    )
    axs[1][:plot](
        ϕ₀,
        -θ₀ .+ 0.9*i,
        linewidth=0.5,
        color="k",
        zorder=-3,
        alpha=0.4
    )
    axs[1].set_xlim(0, 2π)
    m = @. (zoomlim[1] <= ϕ <= zoomlim[2]) | isnan(ϕ)
    axs[2][:plot](
        ϕ[m],
        -θ[m] .+ 0.9*i,
        linewidth=1.5,
        color=colors[i],
        label=rₑ
    )
    axs[2][:plot](
        ϕ₀,
        -θ₀ .+ 0.9*i,
        linewidth=0.5,
        color="k",
        zorder=-3,
        alpha=0.4
    )
    axs[2].set_xlim(zoomlim[1], zoomlim[2])
end
axs[1].set_xticks([0, zoomlim[1], π, zoomlim[2], 2π])
axs[1].set_xticklabels(["0", "", "π", "", "2π"])
axs[2].set_xticks([zoomlim[1], π, zoomlim[2]])
axs[2].set_xticklabels(["15π/16", "π", "17π/16"])
for ax ∈ axs
    ax[:yaxis].set_visible(false)
    ax[:spines][:left].set_visible(false)
    ax[:grid](false)
end
leg = axs[2][:legend]()
leg[:set_title]("Ejecta Multiple")
axs[end][:set_xlabel]("Longitude [rad]")
fig[:tight_layout]()