push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/DiskJockey/")

using LittleMC

using Distributions
using PDMats

const mat = PDiagMat([2.0^2, 3.0^2])
const dist = DiagNormal([0., 0.], mat)

function Gauss(x::Vector{Float64})
    return logpdf(dist, x)
end


using PDMats

jump = full(PDiagMat([1.5^2, 1.5^2]))

println(typeof(jump))

mc = MC(Gauss, 100000, [1.0, 1.0], jump)

start(mc)

println(mean(mc.samples, 2))
println(std(mc.samples, 2))

runstats(mc)

write(mc, "mc.hdf5")

# # Plot the samples
# import PyPlot.plt
# fig, ax = plt.subplots(nrows=2, sharex=true)
# ax[1][:plot](mc.samples'[:, 1])
# ax[2][:plot](mc.samples'[:, 2])
# plt.savefig("plots/littlemc.png")
