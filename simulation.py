import os
import sys
import math
import matplotlib.pyplot as pl
from mpl_toolkits import mplot3d
import numpy as np
import time as tm
import random
from functions import *
from constants import *


ma, mb, mp, radius_a, radius_b, radius_p, a_s, a_p, per_s, per_p, e_s, e_p, i_s, i_p, Omega_s, Omega_p, w_s, w_p = input(
    'kepler16b.txt')

a_a = (mb) / (ma + mb) * a_s  # Calculates semimajor axis of Star A
a_b = (ma) / (ma + mb) * a_s  # Calculates semimajor axis of Star B
n_periods = 3 # number of orbital periods the simulation runs for
f_s, rho_s, res = time(per_s, e_s, per_p, n_periods)  # Import true anomaly and rho from time function
f_p, rho_p, res = time(per_p, e_p, per_p, n_periods)
X_a = [0] * res
Y_a = [0] * res
Z_a = [0] * res
X_b = [0] * res
Y_b = [0] * res
Z_b = [0] * res
X_p = [0] * res
Y_p = [0] * res
Z_p = [0] * res

for i in range(0, res):
    X_a[i] = a_a * rho_s[i] * (math.cos(Omega_s) * math.cos(w_s + f_s[i]) -
                            math.sin(Omega_s) * math.sin(w_s + f_s[i]) * math.cos(i_s))
    Y_a[i] = a_a * rho_s[i] * (math.sin(Omega_s) * math.cos(w_s + f_s[i]) +
                            math.cos(Omega_s) * math.sin(w_s + f_s[i]) * math.cos(i_s))
    Z_a[i] = a_a * rho_s[i] * (math.sin(w_s + f_s[i]) * math.sin(i_s))
    X_b[i] = -a_b * rho_s[i] * (math.cos(Omega_s) * math.cos(w_s + f_s[i]) -
                            math.sin(Omega_s) * math.sin(w_s + f_s[i]) * math.cos(i_s))
    Y_b[i] = -a_b * rho_s[i] * (math.sin(Omega_s) * math.cos(w_s + f_s[i]) +
                            math.cos(Omega_s) * math.sin(w_s + f_s[i]) * math.cos(i_s))
    Z_b[i] = -a_b * rho_s[i] * (math.sin(w_s + f_s[i]) * math.sin(i_s))
    X_p[i] = a_p * rho_p[i] * (math.cos(Omega_p) * math.cos(w_p + f_p[i]) -
                            math.sin(Omega_p) * math.sin(w_p + f_p[i]) * math.cos(i_p))
    Y_p[i] = a_p * rho_p[i] * (math.sin(Omega_p) * math.cos(w_p + f_p[i]) +
                            math.cos(Omega_p) * math.sin(w_p + f_p[i]) * math.cos(i_p))
    Z_p[i] = a_p * rho_p[i] * (math.sin(w_p + f_p[i]) * math.sin(i_p))


d_ab = [0,0]
d_ap = [0,0]
d_bp = [0,0]
signature_a = [0] * n_periods # signature (code) for planet over star A
signature_b = [0] * n_periods # signature (code) for planet over star B
per_p_timestep = res / n_periods # one orbital period in units of time step
for i in range(0,res): 
	d_ap[1] = eucdist(X_a[i], X_p[i], Y_a[i], Y_p[i]) # distance in XY plane of A,P
	d_bp[1] = eucdist(X_b[i], X_p[i], Y_b[i], Y_p[i]) # distance in XY plane of B,P
	for j in range(0, n_periods):
		if (j * per_p_timestep < i and i < (j+1) * per_p_timestep):
			if (d_ap[1] <= radius_a + radius_p and Z_p[i] > Z_a[i] and d_ap[0] > radius_a + radius_p):
				signature_a[j] += 1
			if (d_bp[1] <= radius_b + radius_p and Z_p[i] > Z_b[i] and d_bp[0] > radius_b + radius_p):
				signature_b[j] += 1
	d_ap[0] = d_ap[1]
	d_bp[0] = d_bp[1]
			



print("signatures of A and B")
print(signature_a)
print(signature_b)
print("======================")



if (1 == 0):
	fig = pl.figure()
	ax = pl.axes(projection="3d")
	ax.plot3D([0.], [0.], [0.], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=5, alpha=0.6)
	ax.plot3D(X_a, Y_a, Z_a, color='r')
	ax.plot3D(X_b, Y_b, Z_b, color='g')
	ax.plot3D(X_p, Y_p, Z_p, color='b')
	ax.set_xlabel("X")
	ax.set_ylabel("Y")
	ax.set_zlabel("Z")
	ax.set_xlim(-1.2e11, 1.2e11)
	ax.set_ylim(-1.2e11, 1.2e11)
	ax.set_zlim(-1.2e11, 1.2e11)
	pl.show()

# print("star stats: ecc, inc, Omega, w")
# print(e_s)
# print(i_s)
# print(Omega_s)
# print(w_s)

'''
Notes for tomorrow: 
Today, finished the orbit sims and everything should be running smoothly. 

Tomorrow, create a variable for number of orbital periods the simulation is run for. 

Create the transit detection system, throw out occultations. 
[====DONE ABOVE====]
Then maybe work on the fix all system parameters except for 
the initial phase of the binary and planet. Figure out a way to randomize all 
variables to set up to start running large-scale simulations. 

Create the binary detection system, make some way to store the system parameters
and organize by the detection code. Once there, running simulations and such seems
like the logical next step. 
'''