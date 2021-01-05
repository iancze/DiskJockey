# Download and install RADMC-3D
binary_file = "master.zip"

# Change to the src directory
cd("src")

# Download and extract the RADMC-3D executable
download("https://github.com/dullemond/radmc3d-0.41/archive/$binary_file", binary_file)

# Remove the directory, if a stale version exists from a previous build
if ispath("radmc-3d")
    rm("radmc-3d", recursive=true)
end
# Unzipping the file
run(`unzip $binary_file`)

# make the binary
cd("radmc3d-0.41-master/src")
run(`make`)
