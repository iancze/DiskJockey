# Download and install RADMC-3D
binary_file = "radmc-3d_v0.41_07.07.17.zip"

# Change to the src directory
cd("src")
println("Current working directory ", pwd())

# Download and extract the RADMC-3D executable
download("http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/$binary_file", binary_file)
println("Current working directory ", pwd())

# Unzipping the file
run(`unzip $binary_file`)

cd("radmc-3d/version_0.41/src")
run(`make`)
