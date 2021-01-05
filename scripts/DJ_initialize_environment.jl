#!/usr/bin/env julia

using Pkg

# add the packages to the current environment, defined by JULIA_PROJECT 
Pkg.add(PackageSpec(url = "https://github.com/iancze/DiskJockey.git"))
Pkg.add("ArgParse")
Pkg.add("HDF5")
Pkg.add("PyPlot")
Pkg.add("LaTeXStrings")
Pkg.add("YAML")
Pkg.add("Printf")
Pkg.add("PyCall")
Pkg.add("Images")
Pkg.add("NaNMath")
Pkg.add("NPZ")
Pkg.add("Statistics")

Pkg.instantiate()