import math
from functions import *
from system import *

if (True):
	cbs = System(mu=0.3,e_s=0.3,e_p=0.01,r_a=0.001,r_b=0.003,r_p=0.0004,deltaI=0)
	cbs.simulation()
	cbs.transits()
	cbs.plotOrbits()