using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using DataFrames
using CSV
using PyPlot

pygui(true)

## load the big table of results

df = CSV.read(datadir("sims", "simulations.csv"), DataFrame)

## pick latitude and overlap

θₛ = π/4
Δ = 5.0

##

figure()
for (i,rₑ) ∈ enumerate(unique(df[!,:re]))
    #filter
    dg = df[(df.re .≈ rₑ) .& (df.theta .≈ θₛ) .& (df.overlap .≈ Δ),:]
    #time axis
    t = dg[!,:t]
    #quantiles

    plot(t, dg[!,:f25], color='C'*string(i-1), alpha=0.5, linewidth=0.5, linestyle=":")
    plot(t, dg[!,:f50], color='C'*string(i-1), label="ejecta multiple = $rₑ", linewidth=1.5)
    plot(t, dg[!,:f75], color='C'*string(i-1), alpha=0.5, linewidth=0.5, linestyle=":")
end
xlim(maximum(df[!,:t]), minimum(df[!,:t]))
xlabel("Time [Ga]")
ylabel("Fraction of Shoreline Destroyed")
title("Shoreline Latitude = $(Int(round((π/2 - θₛ)*180/π))) degrees")
legend()
tight_layout()
savefig(plotsdir("fraction_destroyed.png"), dpi=750)

##

figure()
for (i,rₑ) ∈ enumerate(unique(df[!,:re]))
    dg = df[(df.re .≈ rₑ) .& (df.theta .≈ θₛ) .& (df.overlap .≈ Δ),:]
    t = dg[!,:t]
    plot(t, dg[!,:smax25]/1e3, color='C'*string(i-1), alpha=0.5, linewidth=0.5, linestyle=":")
    plot(t, dg[!,:smax50]/1e3, color='C'*string(i-1), label="ejecta multiple = $rₑ", linewidth=1.5)
    plot(t, dg[!,:smax75]/1e3, color='C'*string(i-1), alpha=0.5, linewidth=0.5, linestyle=":")
end
xlim(maximum(df[!,:t]), minimum(df[!,:t]))
xlabel("Time [Ga]")
ylabel("Maximum Shoreline Segment Length [km]")
title("Shoreline Latitude = $(Int(round((π/2 - θₛ)*180/π))) degrees")
legend()
tight_layout()
savefig(plotsdir("maximum_segment_length.png"), dpi=750)

##

figure()
for (i,rₑ) ∈ enumerate(unique(df[!,:re]))
    dg = df[(df.re .≈ rₑ) .& (df.theta .≈ θₛ) .& (df.overlap .≈ Δ),:]
    t = dg[!,:t]
    plot(t, dg[!,:smean25]/1e3, color='C'*string(i-1), alpha=0.5, linewidth=0.5, linestyle=":")
    plot(t, dg[!,:smean50]/1e3, color='C'*string(i-1), label="ejecta multiple = $rₑ", linewidth=1.5)
    plot(t, dg[!,:smean75]/1e3, color='C'*string(i-1), alpha=0.5, linewidth=0.5, linestyle=":")
end
xlim(maximum(df[!,:t]), minimum(df[!,:t]))
xlabel("Time [Ga]")
ylabel("Mean Shoreline Segment Length [km]")
title("Shoreline Latitude = $(Int(round((π/2 - θₛ)*180/π))) degrees")
legend()
tight_layout()
savefig(plotsdir("mean_segment_length.png"), dpi=750)