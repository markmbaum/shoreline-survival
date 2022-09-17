# Shoreline Survival 

## Associated Info  

* A conference abstract describing the study that relies on this repository (as of Jan 5, 2022) can be downloaded from this repository [here](papers/LPSC_2022_abstract.pdf).
* The submitted (March 17, 2022) manuscript can also be downloaded from this repository [here](papers/paper_submitted.pdf)
* The submitted manuscript is also [on arXiv](https://arxiv.org/abs/2206.09816).
* physicsworld magazine did a [short article on the study](https://physicsworld.com/a/simulations-show-that-asteroid-impacts-would-destroy-evidence-for-relic-shorelines-on-mars/).
* The final manuscript was [published by Icarus](https://linkinghub.elsevier.com/retrieve/pii/S0019103522002792).
* This code is also archived at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6354413.svg)](https://doi.org/10.5281/zenodo.6354413)
* Data files associated with this project are archived at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6354426.svg)](https://doi.org/10.5281/zenodo.6354426)

## Code

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make a reproducible scientific project named

> Shoreline Survival

It is authored by Mark Baum <markmbaum@protonmail.com>.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the

git-history and may need to be downloaded independently.

1. Open a Julia console and do:

```

julia> using Pkg
julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
julia> Pkg.activate("path/to/this/project")
julia> Pkg.instantiate()

```

This will install all necessary packages for you to be able to run the scripts and everything should work out of the box, including correctly finding local paths.
