import os
import sys
import math
import matplotlib.pyplot as pl
from mpl_toolkits import mplot3d
import random
import numpy as np
from functions import *

class System:
	def __init__(self,mu,e_s,e_p,r_a,r_b,r_p,deltaI,n_per=3,res=10**5):
		self.mu = mu
		self.e_s = e_s
		self.e_p = e_p
		self.r_a = r_a
		self.r_b = r_b
		self.r_p = r_p
		self.i_s = math.pi/2
		self.i_p = self.i_s + deltaI
		self.M0 = random.uniform(0,2*math.pi)
		self.Omega_s = 0
		self.Omega_p = random.uniform(-0.01*math.pi, 0.01*math.pi)
		self.w_s = random.uniform(0,2*math.pi)
		self.w_p = random.uniform(0,2*math.pi)
		self.a_a = self.mu
		self.a_b = 1-self.mu
		self.m_a = self.mu
		self.m_b = 1-self.mu
		self.n_per = n_per
		self.res = res

		return None

	def time(self,e,per):
		'''
		time: Given binary eccentricity and period of the body, computes true anomaly and rho.
		rho = (1-e^2)/(1+e*cos(f)), r=a*rho

		e: binary eccentricity
		per_p: orbital period, used for time scale
		n_per: number of periods, used for time scale. Default: 3
		res: resolution (number of points) Default: 1e5

		t_range: time range. 
		n: mean motion (2pi/P)
		dt: time step
		M: mean anomaly
		E: eccentric anomaly
		Eold, deltaE: dummy variables for iteration
		'''
		t_range = self.n_per * per
		tol = 1e-6
		n = 2*math.pi / per
		dt = n * t_range / self.res
		M = np.linspace(self.M0,self.M0+n*t_range,self.res)
		E = [3.14] * self.res
		Eold = 0
		for i in range(self.res):
			deltaE = 1
			while (deltaE>tol):
				Eold = E[i]
				E[i] = E[i] - (E[i] - e * math.sin(E[i]) - M[i]) / (1 - e * math.cos(E[i]))
				deltaE = abs(E[i]-Eold)
		f = [0] * self.res
		rho = [0] * self.res
		for i in range(self.res):
			f[i] = 2 * math.atan(((1 + e) / (1 - e))**0.5 * math.tan(E[i] / 2))
		for i in range(self.res):
			rho[i] = (1 - e**2) / (1 + e * math.cos(f[i]))
		return f, rho

	def simulation(self):
		'''
		simulation: Runs the simulation with given parameters.

		a_p: semimajor axis of planet, defined as 1.2 times the HW99 critical semimajor axis.
		per_p: period of the planet.
		f: true anomaly (of star or planet)
		rho: normalized distance of body (r = a * rho)
		res: resolution; number of points. Default: 10^5.
		X,Y,Z: reference coordinates. (of star A, star B, and planet P)
		'''
		self.a_p = 1.2*HW99(self.e_s,self.mu)
		self.per_p = self.a_p**1.5
		self.f_s,self.rho_s = self.time(self.e_s,1)
		self.f_p,self.rho_p = self.time(self.e_p,self.per_p)
		self.X_a = [0] * self.res
		self.Y_a = [0] * self.res
		self.Z_a = [0] * self.res
		self.X_b = [0] * self.res
		self.Y_b = [0] * self.res
		self.Z_b = [0] * self.res
		self.X_p = [0] * self.res
		self.Y_p = [0] * self.res
		self.Z_p = [0] * self.res
		for i in range(self.res):
			self.X_a[i] = self.a_a * self.rho_s[i] * (math.cos(self.Omega_s) * math.cos(self.w_s + self.f_s[i]) -
									math.sin(self.Omega_s) * math.sin(self.w_s + self.f_s[i]) * math.cos(self.i_s))
			self.Y_a[i] = self.a_a * self.rho_s[i] * (math.sin(self.Omega_s) * math.cos(self.w_s + self.f_s[i]) +
									math.cos(self.Omega_s) * math.sin(self.w_s + self.f_s[i]) * math.cos(self.i_s))
			self.Z_a[i] = self.a_a * self.rho_s[i] * (math.sin(self.w_s + self.f_s[i]) * math.sin(self.i_s))
			self.X_b[i] = -self.a_b * self.rho_s[i] * (math.cos(self.Omega_s) * math.cos(self.w_s + self.f_s[i]) -
									math.sin(self.Omega_s) * math.sin(self.w_s + self.f_s[i]) * math.cos(self.i_s))
			self.Y_b[i] = -self.a_b * self.rho_s[i] * (math.sin(self.Omega_s) * math.cos(self.w_s + self.f_s[i]) +
									math.cos(self.Omega_s) * math.sin(self.w_s + self.f_s[i]) * math.cos(self.i_s))
			self.Z_b[i] = -self.a_b * self.rho_s[i] * (math.sin(self.w_s + self.f_s[i]) * math.sin(self.i_s))
			self.X_p[i] = self.a_p * self.rho_p[i] * (math.cos(self.Omega_p) * math.cos(self.w_p + self.f_p[i]) -
									math.sin(self.Omega_p) * math.sin(self.w_p + self.f_p[i]) * math.cos(self.i_p))
			self.Y_p[i] = self.a_p * self.rho_p[i] * (math.sin(self.Omega_p) * math.cos(self.w_p + self.f_p[i]) +
									math.cos(self.Omega_p) * math.sin(self.w_p + self.f_p[i]) * math.cos(self.i_p))
			self.Z_p[i] = self.a_p * self.rho_p[i] * (math.sin(self.w_p + self.f_p[i]) * math.sin(self.i_p))
		return None
	def transits(self):
		'''
		transits: Calculates the transits per period. Prints the planet/star A and planet/star B transits.

		signature: array with n_per entries, containing the number of transits in each period
		currently_transiting: boolean, True if planet is currently transiting star A/B.
		current_period: current period, from 0 to n_per-1.
		xyz: array ([x,y,z]) containing xyz coordinates. 
		points_per_period: number of data points per orbital period. 
		'''
		self.signature_a = [0]*self.n_per
		self.signature_b = [0]*self.n_per
		currently_transitingA = False
		currently_transitingB = False
		current_period = 0
		xyza = [0,0,0]
		xyzb = [0,0,0]
		xyzp = [0,0,0]
		points_per_period = int(self.res/self.n_per)
		for i in range(self.res):
			xyza = [self.X_a[i],self.Y_a[i],self.Z_a[i]]
			xyzb = [self.X_b[i],self.Y_b[i],self.Z_b[i]]
			xyzp = [self.X_p[i],self.Y_p[i],self.Z_p[i]]
			if (i % points_per_period == points_per_period-1):
				current_period += 1
			if (currently_transitingA and not isTransiting(xyza,xyzp,self.r_a,self.r_p)):
				currently_transitingA = False
			elif (not currently_transitingA and isTransiting(xyza,xyzp,self.r_a,self.r_p)):
				currently_transitingA = True
				self.signature_a[current_period] += 1
			if (currently_transitingB and not isTransiting(xyzb,xyzp,self.r_b,self.r_p)):
				currently_transitingB = False
			elif (not currently_transitingB and isTransiting(xyzb,xyzp,self.r_b,self.r_p)):
				currently_transitingB = True
				self.signature_b[current_period] += 1
		print("Planet-Star A transits: " + str(self.signature_a))
		print("Planet-Star B transits: " + str(self.signature_b))
		return None

	def plotOrbits(self, center_of_mass=True, scale=2.5, colorA='r', colorB='g', colorP='b'):
		'''
		plotOrbits: Plots 3D orbits of star A, star B, and planet P. 

		center_of_mass: boolean. If true, plots a point to display the center of mass. Default: True
		scale: sets the scale of the axes. Default: 2.5
		color: colors of each path. Default: red A, green B, blue P.
		'''
		fig = pl.figure()
		ax = pl.axes(projection="3d")
		if(center_of_mass):
			ax.plot3D([0.], [0.], [0.], markerfacecolor='k', markeredgecolor='k', marker='o', markersize=5, alpha=0.6)
		ax.plot3D(self.X_a, self.Y_a, self.Z_a, color=colorA)
		ax.plot3D(self.X_b, self.Y_b, self.Z_b, color=colorB)
		ax.plot3D(self.X_p, self.Y_p, self.Z_p, color=colorP)
		ax.set_xlabel("X")
		ax.set_ylabel("Y")
		ax.set_zlabel("Z")
		ax.set_xlim(-scale,scale)
		ax.set_ylim(-scale,scale)
		ax.set_zlim(-scale,scale)
		ax.azim=-90
		ax.elev=90
		pl.show()
		return None