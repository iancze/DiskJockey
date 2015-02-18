# Go through and add all required packages

pkgs = ["YAML", "PyPlot", "LaTeXStrings", "Distributions", "PDMats", "HDF5", "ArgParse", "NPZ", "Logging"]

for pkg in pkgs
    Pkg.add(pkg)
end
