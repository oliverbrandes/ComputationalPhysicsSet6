import numpy as np
import scipy
import scipy.linalg as lag

class EllipticSolver:

    def __init__(self, r):
        self.r = r
        self.N = len(r)
        self.rhs = np.zeros(self.N)
        self.AB = np.zeros((3, self.N))
        self.delta = r[1]-r[0]


    def setup_matrix(self):
        N, r = self.N, self.r
        delta = self.delta

         #diagonals
        l_diag = 1/(delta**2) - 1/(r*delta)
        diag = -2/(delta**2) * np.ones(N)
        t_diag = 1/(delta**2) + 1/(r*delta)

         #lower boundary as specified by problem 4.3b
        diag[0] = -3/(delta**2)  #A[0][0]
        t_diag[0] = 3/(delta**2) #A[0][1]

         #upper boundary as specified by 4.3c
        diag[-1] = -1/(delta**2)
        l_diag[-1] = (N - 3/2)/((N-1/2)*delta**2)

         #storing the diagonals in matrix AB
        self.AB[0, 1:] = t_diag[:-1]
        self.AB[1, :] = diag
        self.AB[2, :-1] = l_diag[1:]



    def setup_rhs(self,factor,fct):
        self.rhs = factor*fct

    def solve(self):
        sol = lag.solve_banded((1,1), self.AB, self.rhs)
        return sol