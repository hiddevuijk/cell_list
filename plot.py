import numpy as np
import matplotlib.pyplot as plt
from sys import exit

Lx = 10.0
Ly = 10.0
positions = np.loadtxt("positions.dat")

for pi in positions:
    x = pi[0]
    y = pi[1]
    x -= Lx * np.floor(x / Lx)
    y -= Ly * np.floor(y / Ly)
    plt.scatter(x,y)

plt.show()

