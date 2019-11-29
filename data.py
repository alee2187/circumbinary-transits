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
'''
Parameters needed: 
ma, mb, mp, radius_a, radius_b, radius_p, a_s, a_p, per_s, per_p, 
e_s, e_p, i_s, i_p, Omega_s, Omega_p, w_s, w_p, M0_s, M0_p

Omega_s = 0
Parameters to generate: 
mu, radius_a, radius_b, radius_p, a_s, e_s, e_p, i_s, i_p, Omega_p, w_s, w_p, M0_s, M0_p

Other parameters to specify: 
n_periods, res, tol



mu: [0,0.5]
radius_a (Rsun): [,]
radius_b (Rsun): [,]
radius_p (Rjup): [0.25,1.05]
a_s (AU): [0.0836,0.22882]
e_s: [0.023,0.521]
e_p: [0.007,0.411]
i_s (deg): [87,93]
i_p (deg): [87,93]
Omega_p (deg): [-2,2]
w_s (deg): [,]
w_p (deg): [,]
M0_s (deg): [0,360]
M0_p (deg): [0,360]
'''
def datain(grid=10):
	mass_a = np.linspace(0.6897, 1.47, grid) # random
	mass_b = np.linspace(0.1951, 1.0208, grid) # random
	radius_a = np.linspace(0.6489, 1.79, grid) # random
	radius_b = np.linspace(0.2143, 1.0927, grid) # random
	radius_p = np.linspace(0.25, 1.05, grid) # random
	a_s = np.linspace(0.0836, 0.22882, grid) # random
	e_s = np.linspace(0.023,0.521,grid) # set
	e_p = np.linspace(0.007,0.411,grid) # set
	i_s = np.linspace(87,93,grid) # set
	i_p = np.linspace(87,93,grid) # set
	Omega_p = np.linspace(-2,2,grid) # random
	w_s = np.linspace(0,359,grid) # random
	w_p = np.linspace(0,359,grid) # random
	M0_s = np.linspace(0,359,grid) # Random
	M0_p = np.linspace(0,359,grid) # Random

	mass_a *= MSUN_KG
	mass_b *= MSUN_KG
	radius_a *= RSUN_M
	radius_b *= RSUN_M
	radius_p *= RJUPITER_M
	a_s *= AU_M
	i_s *= math.pi/180
	i_p *= math.pi/180
	Omega_p *= math.pi/180
	w_s *= math.pi/180
	w_p *= math.pi/180
	M0_s *= math.pi/180
	M0_p *= math.pi/180
	data_out = (mass_a, mass_b, radius_a, radius_b, radius_p, a_s, e_s, e_p, i_s, i_p, Omega_p, w_s, w_p, M0_s, M0_p)
	return data_out

grid = 10
data = datain(grid)
print(data)
x=0
for i in range (0,10):
	print(x)