import scipy
from scipy.interpolate import barycentric_interpolate
import numpy as np
import matplotlib.pyplot as plt

#create an array between -1 and 1 with N points, this will actually create an array


def even_spaced_array(N):
    delta = 2.0/(N-1)
    sol_array = np.zeros(N)
    for i in range(N):
        sol_array[i] = -1 + delta * i
    return sol_array

#create user defined function (function not attached to class (apparently)) for
# fucntion defined in 6_1

def fct(x:float):
    return np.exp(-x**2)/(1+25*x**2)

#creating array based on x for y using f(x) (as seen in 6_1)

def y_array(x):
    sol_array = np.zeros(len(x))
    for i, x in enumerate(x):
        sol_array[i] = fct(x)
    return sol_array

#creating bary_centric interpolation function on already created array

def bary_interpolates(x, y, x_int):
    sol = barycentric_interpolate(x, y, x_int)
    return sol

#creating x values of a centric grid for N = 5 and N = 20
N_is_5 = even_spaced_array(5)
N_is_20 = even_spaced_array(20)

#this is creating N spaced array that is significantly larger than either N = 5
# or N = 20
x_int = even_spaced_array(200)


plt.plot( x_int, bary_interpolates(N_is_5, y_array(N_is_5), x_int), "-", label="Bary Centric Interpolation for N = 5")
plt.plot(x_int, bary_interpolates(N_is_20, y_array(N_is_20), x_int), "-", label="Bary Centric Interpolation for N = 20")
plt.plot(x_int, fct(x_int), "-", label="Analytical Solution for f(x)")
plt.plot(N_is_5, y_array(N_is_5), "o", label="N = 5")
plt.plot(N_is_20, y_array(N_is_20), "o", label="N = 20")
plt.xlabel("X, evenly spaced intervals between -1 and 1")
plt.ylabel("Y, the interpolated value of function f(x) in 6_1")
plt.title("Bary Centric Interpolation")
plt.legend()
plt.show()



#checks for functions

#check for even_spaced_grid
#even_spaced_array = even_spaced_array(20)
#print(even_spaced_array)


#check for y_array
#y_array = y_array(even_spaced_array)
#print(y_array)

#check for bary_interpolates
#print(bary_interpolates(even_spaced_array))

