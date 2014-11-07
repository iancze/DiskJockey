#Test plotting

import PyPlot.plt

rands = rand((4, 4, 4))

fig, ax = plt.subplots(ncols=4, figsize=(7, 4))

println(ax)

ax[1]

ax[1][:imshow](rands[:,:,1])


plt.savefig("test.png")


