import os
import sys
import math

def HW99(e,mu):
	'''
	HW99: Computes Holman-Wiegert 1999 critical semimajor axis. Returns the ratio
	a_c/a_bin.

	e: binary eccentricity
	mu: mass ratio, in range [0,0.5]
	'''
	return (1.60) + (5.10 * e) + (-2.22 * e**2) + (4.12 * mu) + (-5.09 * mu**2) + (-4.27 * e * mu) + (4.61 * (e * mu)**2)

def dist(x1,x2,y1,y2):
	'''
	dist: computes Euclidean distance.
	'''
	return ((x2-x1)**2 + (y2-y1)**2)**0.5

def isTransiting(xyzs,xyzp,r_s,r_p): 
	'''
	isTransiting: checks if two bodies are transiting and if 'p' is in front of 's'. 
	The observer is placed at Z = +infinity. 
	
	xyzs: array, xyz coordinates of star ([x,y,z])
	xyzp: array, xyz coordinates of planet ([x,y,z])
	r_s: radius of star
	r_p: radius of planet

	xydist: distance between bodies on the XY plane. 
	zdist: distance between bodies on the Z axis. If zdist is positive, 'p' is in front of 's'. 
	covered: boolean, true if bodies are overlapping on the XY plane.
	planet_front: boolean, true if 'p' is in front of 's'.
	'''

	xydist = dist(xyzs[0],xyzp[0],xyzs[1],xyzp[1])
	zdist = xyzp[2]-xyzs[2]
	covered = False
	planet_front = False
	if (xydist < r_s + r_p):
		covered = True
	if (zdist > 0):
		planet_front = True

	return (covered and planet_front)
