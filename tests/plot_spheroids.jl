push!(LOAD_PATH, "/home/ian/Grad/Research/Disks/JudithExcalibur/")

using gridding
import PyPlot.plt
using LaTeXStrings

etas = linspace(-1.0, 1.0, 100)

alphas = Float64[0., 0.5, 1.0, 1.5, 2.0]

fig = plt.figure(figsize=(6,6))
ax = fig[:add_subplot](111)

for alpha in alphas
    label = @sprintf("%.1f", alpha)
    ax[:plot](etas, spheroid(etas, alpha), label=label)
end

ax[:set_xlabel](L"$\eta$")
ax[:legend](loc="upper right")
plt.savefig("../plots/spheroids.png")

fig = plt.figure(figsize=(6,6))
ax = fig[:add_subplot](111)

for alpha in alphas
    label = @sprintf("%.1f", alpha)
    ax[:plot](etas, corrfun(etas, alpha), label=label)
end

ax[:set_xlabel](L"$\eta$")
ax[:legend](loc="upper right")
plt.savefig("../plots/corrfuns.png")

fig = plt.figure(figsize=(6,6))
ax = fig[:add_subplot](111)

for alpha in alphas
    label = @sprintf("%.1f", alpha)
    ax[:plot](etas, gcffun(etas, alpha), label=label)
end

ax[:set_xlabel](L"$\eta$")
ax[:legend](loc="upper right")
plt.savefig("../plots/gcffuns.png")
