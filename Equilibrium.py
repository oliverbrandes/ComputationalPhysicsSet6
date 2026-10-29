
import numpy as np
import scipy


class Hydrostatic_Equilibrium:
    kay: float
    gamma: float

    def __init__(self, kay: float, gamma: float):
        self.kay = kay
        self.gamma = gamma

    def P_of_rho(self, rho):
        return self.kay * rho ** self.gamma

    def rho_of_P(self, P:float):
        rho = (P/self.kay)**(1/self.gamma)
        return rho

    def derivs(self, r, vars):
        P = vars[0]
        m = vars[1]
        rho = self.rho_of_P(P)
        if(r > 0.0):
            dPdr = -(rho)*1*m/(r**2) #G is 1 here
            dmdr = 4*np.pi*rho*r**2
            # dvars = np.array([dPdr, dmdr])
        else:
            dPdr = 0
            dmdr = 0
        dvars = np.array([dPdr, dmdr])
        return dvars

    def reached_surface(self, r, vars):
        return vars[0]

    def equilibrium(self, rho_c):
        Hydrostatic_Equilibrium.reached_surface.terminal = True
        Hydrostatic_Equilibrium.reached_surface.direction = -1
        sol = scipy.integrate.solve_ivp(self.derivs, [0, 10],[self.P_of_rho(rho_c), 0], events = self.reached_surface)
        P_array = sol.y[0]
        m_array = sol.y[1]
        M = m_array[-1]
        rho_array = self.rho_of_P(P_array)
        R = sol.t[-1]
        r_array = sol.t
        return M, R, rho_array, P_array, m_array, r_array

# star = Hydrostatic_Equilibrium(1, 2)   # K=1, Gamma=2
# M, R, rho_array, P_array, m_array, r_array = star.equilibrium(1)   # rho_c = 1
# print(P_array[-1])