using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using DataFrames
using CSV
using PyPlot

pygui(true)

## segments

dirseg = datadir("sims", "mapped-segments")
segments = Dict{Int64,Vector{SphericalSegment}}()
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
for s ∈ segments[idx]
    plot(
        (s.a.ϕ, s.b.ϕ),
        (s.a.θ, s.b.θ),
        color="C0",
        alpha=0.8,
        linewidth=1.5
    )
end