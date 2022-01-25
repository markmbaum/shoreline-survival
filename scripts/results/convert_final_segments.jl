using DrWatson
@quickactivate "Shoreline Survival"
push!(LOAD_PATH, srcdir())
using ShorelineSurvival
using Base.Threads: @threads

##

function segs2txt(path, S)::Nothing
    open(path, "w") do io
        println(io, "theta_a, phi_a, theta_b, phi_b")
        for s ∈ S
            println(io, s.a.θ, ',', s.a.ϕ, ',', s.b.θ, ',', s.b.ϕ)
        end
    end
    return nothing
end

##

dirseg = datadir("sims", "mapped-segments")
@threads for path ∈ readdir(dirseg, join=true)
    if path[end-3:end] != ".txt"
        #read segments
        S = loadsegments(path)
        #write them to a text file with the same name
        segs2txt(path*".txt", S)
    end
end