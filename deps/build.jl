# Download and install RADMC-3D
binary_file = "radmc-3d_v0.41_07.07.17.zip"

# Change to the src directory
cd("src")

# Download and extract the RADMC-3D executable
download("http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/$binary_file", binary_file)

# Remove the directory, if a stale version exists from a previous build
if ispath("radmc-3d")
    rm("radmc-3d", recursive=true)
end
# Unzipping the file
run(`unzip $binary_file`)

# make the binary
cd("radmc-3d/version_0.41/src")
run(`make`)
