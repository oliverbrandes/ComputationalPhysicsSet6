import numpy as np
import matplotlib.pyplot as plt
from InData import InData

K = 1.0
rho_c = 1.0
N = 50

#Gamma = 4/3
r_furthest = 4.5
delta = r_furthest/N
r = (np.arange(N) + 0.5)*delta          # cell-centered grid
M, R, rho, P, m = InData(K, 4/3).equilibrium(rho_c, r)

#checks for M, r, rho, P, m
# print (M, m)
# print (R, rho, P)
# print(sum(m))

plt.plot(r, rho, "o", label="Numerical density on cell centered grid")
plt.xlabel("r, cell centered radial grid points")
plt.ylabel("rho, the density of the star")
plt.title("Gamma = 4/3, rho_c = 1")
plt.legend()
plt.show()

#Gamma = 2
r_furthest = 1.5
delta = r_furthest/N
r = (np.arange(N) + 0.5)*delta          # cell-centered grid

M, R, rho, P, m = InData(K, 2.0).equilibrium(rho_c, r)

#fine grid for the analytical solution (inside the star only)
r_analytical = np.linspace(1e-6, R, 300)      # start just above 0 to avoid 0/0

plt.plot(r, rho, "o", label="Numerical density on cell centered grid")
plt.plot(r_analytical, rho_c*np.sin(np.pi*r_analytical/R)/(np.pi*r_analytical/R), "-", label="Analytical Solution for rho(r)")
plt.xlabel("r, cell centered radial grid points")
plt.ylabel("rho, density of star")
plt.title("Gamma = 2, rho_c = 1")
plt.legend()
plt.show()