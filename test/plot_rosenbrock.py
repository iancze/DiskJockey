import numpy as np
import matplotlib.pyplot as plt

def rosenbrock(x,y):
    a = 1.
    b = 100.
    return -((a - x)**2 + b * (y - x**2)**2)

N = 100
xs = np.linspace(-3, 3, num=N)
ys = np.linspace(-1, 3, num=N)
XX,YY = np.meshgrid(xs, ys)
ZZ = rosenbrock(XX,YY)
mm = np.max(ZZ)
plt.contour(XX,YY, ZZ, levels=np.linspace(mm - 10, mm, num=10))
plt.savefig("contour.png")
