import numpy as np
import matplotlib.pyplot as plt
from InData import InData
from EllipticSolver import EllipticSolver

K = 1.0
rho_c = 1.0
N = 100

for gamma, r_out in [(4/3, 7.0), (2.0, 2.5)]:
    #cell-centered grid
    delta = r_out/N
    r = (np.arange(N) + 0.5)*delta

    #equilibrium model from 6.2
    M, R, rho, P, m = InData(K, gamma).equilibrium(rho_c, r)

    #solve laplacian(phi) = 4 pi rho
    density_distribution = EllipticSolver(r)
    density_distribution.setup_matrix()
    density_distribution.setup_rhs(4*np.pi, rho)
    phi = density_distribution.solve()

    #part (c): central value from part (b)
    phi_0 = 9/8*phi[0] - 1/8*phi[1]
    phi_0_alternative = -(gamma/(gamma-1))-M/R
    print(phi_0, phi_0_alternative)
    # print(f"Gamma = {gamma:.3f}:  phi(0) = {phi_0:.4f}")

    plt.plot(r, phi, "o", label="r$\phi$ on cell centered grid")
    plt.plot(0, phi_0, "*", markersize=12, label=f"phi(0) = {phi_0:.4f}")
    plt.xlabel("r, cell centered radial grid points")
    plt.ylabel(r"$\phi$, the gravitational potential")
    plt.title(r"$\phi$ vs r")
    plt.legend()
    plt.show()

