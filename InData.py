import numpy as np
import scipy
import scipy.integrate
import scipy.interpolate


class InData:
    kay: float
    gamma: float

    def __init__(self, kay: float, gamma: float):
        self.kay = kay
        self.gamma = gamma

    def P_of_rho(self, rho):
        return self.kay * rho ** self.gamma

    def rho_of_P(self, P: float):
        P = np.maximum(P, 0)   # P < 0 just past the surface would give nan
        rho = (P/self.kay)**(1/self.gamma)
        return rho

    def derivs(self, r, vars):
        P = vars[0]
        m = vars[1]
        rho = self.rho_of_P(P)
        if(r > 0.0):
            dPdr = -(rho)*1*m/(r**2) #G is 1 here
            dmdr = 4*np.pi*rho*r**2
        else:
            dPdr = 0
            dmdr = 0
        dvars = np.array([dPdr, dmdr])
        return dvars

    def reached_surface(self, r, vars):
        return vars[0]
    reached_surface.terminal = True
    reached_surface.direction = -1


    #simplifying later interpolation step with function
    def interpolate(self, x, y, r, outside_value):
        f = scipy.interpolate.interp1d(x, y, "cubic", bounds_error=False, fill_value = outside_value)
        return f(r)

    def equilibrium(self, rho_c, r):
        delta = r[1] - r[0]
        r_out = r[-1] + delta/2

        #First instance of solve_ivp used terminate at the surface to get R and M
        sol = scipy.integrate.solve_ivp(self.derivs, [0, 10*r_out],[self.P_of_rho(rho_c), 0], events=self.reached_surface)
        R = sol.t[-1]
        M = sol.y[1][-1]

        #Check that R < r_out
        if R >= r_out:
            raise ValueError(f"R = {R:.4f} is not smaller than r_out = {r_out:.4f}")

        #Second instance of solve_ivp evaluate the solution at t_eval (many more points than r)
        t_eval = np.linspace(0, R, 20*len(r))
        sol = scipy.integrate.solve_ivp(self.derivs, [0, R], [self.P_of_rho(rho_c), 0], t_eval=t_eval)
        P_inter = sol.y[0]
        m_inter = sol.y[1]
        rho_inter = self.rho_of_P(P_inter)

        #Interpolate to r, with rho = P = 0 for r > R
        rho_function = self.interpolate(sol.t, rho_inter, r,0)
        P_function = self.interpolate(sol.t, P_inter, r, 0)
        m_function = self.interpolate(sol.t, m_inter, r, 0)

        return M, R, rho_function, P_function, m_function