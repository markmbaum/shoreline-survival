# Shoreline Survival [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5821983.svg)](https://doi.org/10.5281/zenodo.5821983)


This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
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

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

-----

## More Information

* A conference abstract describing the study that relies on this repository (as of Jan 5, 2022) can be downloaded here: [papers/LPSC_2022_abstract.pdf](papers/LPSC_2022_abstract.pdf)
* Data files associated with this project can be found here: [10.5281/zenodo.5821983](https://doi.org/10.5281/zenodo.5821983)
