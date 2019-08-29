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
from data import *

loopout = open('looper.txt', 'w')

grid = 10
dat = datain(grid)
print(dat[0][0])
looper = []
for a1 in range(0,grid):
	for a2 in range(0,grid):
		for a3 in range(0,grid):
			for a4 in range(0,grid):
				for a5 in range(0,grid): 
					for a6 in range(0,grid):
						for a7 in range(0,grid):
							for a8 in range(0,grid):
								for a9 in range(0,grid):
									print(str(10*a1 + a2 + 0.1*a3 + 0.01*a4 + 0.001*a5 + 0.0001*a6 + 0.00001*a7 + 0.000001*a8 + 0.0000001*a9) + "%")
									for a10 in range(0,grid):
										for a11 in range(0,grid):
											for a12 in range(0,grid):
												for a13 in range(0,grid):
													for a14 in range(0,grid):
														for a15 in range(0,grid):
															loopout.write(str(a1) + ' ' +str(a2) + ' ' + str(a3) + ' ' + str(a4) + ' ' + str(a5) + ' ' + str(a6) + ' ' + str(a7) + ' ' + str(a8) + ' ' + str(a9) + ' ' + str(a10) + ' ' + str(a11) + ' ' + str(a12) + ' ' + str(a13) + ' ' + str(a14) + ' ' + str(a15) + '\n')

print(looper)
# sa, sb = simulation(mass_a, mass_b, radius_a, radius_b, radius_p, 
# 	a_s, e_s, e_p, i_s, i_p, Omega_p, w_s, w_p, M0_s, M0_p, n_periods=10)
# print(sa)
# print(sb)
# print()