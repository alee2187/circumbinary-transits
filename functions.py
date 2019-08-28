from __future__ import division, print_function, absolute_import, unicode_literals
import os
import sys
import math
from six import string_types
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# import planetplanet as pp
# from planetplanet import Planet, Star, Moon, System
# from planetplanet.constants import *
# from planetplanet.photo.maps import UniformMap
# from planetplanet import maps
import matplotlib.pyplot as pl
import matplotlib
import numpy as np
import time as tm
import random 
from constants import *
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection


def G(m, l, t):  # Value of Universal Gravitational Constant G
    G_const = 6.67408 * 10**-11
    if (m == 'kg' and l == 'm' and t == 's'):
        G_const = 6.67408 * 10**-11
    elif (m == 'solar' and l == 'solar' and t == 'days'):
        G_const = 2942
    elif (m == 'solar' and l == 'm' and t == 'days'):
        G_const = 9.91 * 10**29
    elif (m == 'solar' and l == 'km' and t == 'days'):
        G_const = 9.91 * 10**20
    elif (m == 'solar' and l == 'AU' and t == 'days'):
        G_const = 2.959 * 10**-4
    elif (m == 'earth' and l == 'earth' and t == 'days'):
        G_const = 11470
    elif (m == 'earth' and l == 'm' and t == 'days'):
        G_const = 2.975 * 10**24
    elif (m == 'earth' and l == 'km' and t == 'days'):
        G_const = 2.975 * 10**15
    elif (m == 'earth' and l == 'AU' and t == 'days'):
        G_const = 8.887 * 10**-10
    return G_const


def a_to_P(a, m, G):  # Semimajor Axis to Period, Kepler's 3rd
    P = (2 * math.pi) / ((m * G)**0.5) * a**1.5
    return P


def P_to_a(P, m, G):  # Period to Semimajor Axis, Kepler's 3rd
    a = ((G * m) / (4 * (math.pi**2)))**(1 / 3) * (P**(2 / 3))
    return a

def HW99(e, mu):
    '''
    HW99: Computes Holman-Wiegert 1999 critical semimajor axis. Returns ratio a_c / a_bin.

    e: binary eccentricity
    mu: mass ratio, in range [0,0.5]

    r: ratio a_c/a_bin
    '''
    r = (1.60) + (5.10 * e) + (-2.22 * e**2) + (4.12 * mu) + \
        (-5.09 * mu**2) + (-4.27 * e * mu) + (4.61 * (e * mu)**2)
    return r

def mu(mass_a, mass_b):
    mu = mass_a / (mass_a + mass_b)
    mu = min(mu, 1 - mu)
    return mu

def input(data=""):
    '''
    input: Provides space for entering parameters for a stellar body.
    Unless specified otherwise, units are MKS. Orbital parameters are in degrees.

    data: asks to import a txt file with the data for a certain system.

    mass_a: mass of star A
    mass_b: mass of star B
    mass_p: mass of planet P
    radius_a: radius of star A
    radius_b: radius of star B
    radius_p: radius of planet P
    a_s: binary separation; a_a + a_b
    a_p: semimajor axis of planet
    per_s: period of binary stars
    per_p: period of planet
    e_s: binary eccentricity
    e_p: orbital eccentricity of the planet
    i_s: inclination of the stars
    i_p: inclination of the planet
    Omega_s: Longitude of the ascending node of the stars
    Omega_p: Longitude of the ascending node of the planets
    w_s: Argument of periapsis of the stars
    w_p: Argument of periapsis of the planet

    '''
    dummy = 0
    if (data == ""):
        dummy = 1
        return 0
    else:
        fin = open(data, 'r')

        mass_a = float(fin.readline()) * MSUN_KG
        mass_b = float(fin.readline()) * MSUN_KG
        mass_p = float(fin.readline()) * MJUPITER_KG
        radius_a = float(fin.readline()) * RSUN_M
        radius_b = float(fin.readline()) * RSUN_M
        radius_p = float(fin.readline()) * RJUPITER_M
        a_s = float(fin.readline()) * AU_M
        a_p = float(fin.readline()) * AU_M
        per_s = float(fin.readline()) * DAY_S
        per_p = float(fin.readline()) * DAY_S
        e_s = float(fin.readline())
        e_p = float(fin.readline())
        i_s = float(fin.readline()) * RAD_DEG
        i_p = float(fin.readline()) * RAD_DEG
        Omega_s = float(fin.readline()) * RAD_DEG
        Omega_p = float(fin.readline()) * RAD_DEG
        w_s = float(fin.readline()) * RAD_DEG
        w_p = float(fin.readline()) * RAD_DEG
        return mass_a, mass_b, mass_p, radius_a, radius_b, radius_p, a_s, a_p, per_s, per_p, e_s, e_p, i_s, i_p, Omega_s, Omega_p, w_s, w_p


def time(per, e, orbper, n_periods):
    '''
    time: Given period and eccentricity, along with the orbital period of the planet, 
    computes true anomaly, rho, and outputs resolution.
    rho = (1 - e^2) / (1 + e * cos f), r = a * rho

    per: period of the body
    e: eccentricity of the body
    orbper: orbital period of the planet, used for time scale
    n_periods: number of periods, used for the time scale

    t_range: time range
    res: resolution (number of data points)
    M0: initial mean anomaly
    tol: tolerance for Newton interation when solving for eccentric anomaly
    n: variable defined for mean anomaly equation
    dt: time step
    M: mean anomaly 
    E: eccentric anomaly
    Eold: dummy variable for iteration
    deltaE: dummy variable for iteration

    '''
    t_range = n_periods * orbper
    res = 10**5
    M0 = 0
    tol = 1e-6
    n = 2 * math.pi / per
    dt = n * t_range / res
    M = np.linspace(M0, M0 + n * t_range, res)
    E = [3.14] * res
    Eold = 0
    for i in range(0, res):
        deltaE = 1
        while (deltaE > tol):
            Eold = E[i]
            E[i] = E[i] - (E[i] - e * math.sin(E[i]) - M[i]) / \
                (1 - e * math.cos(E[i]))
            deltaE = abs(E[i] - Eold)
            # print(deltaE)
        # print(i)
    f = [0] * res
    rho = [0] * res
    for i in range(0, res):
        f[i] = 2 * math.atan(((1 + e) / (1 - e))**0.5 * math.tan(E[i] / 2))
    for i in range(0, res):
        rho[i] = (1 - e**2) / (1 + e * math.cos(f[i]))
    return f, rho, res

def eucdist(x1, x2, y1, y2):
    '''
    eucdist: computes Euclidean distance.
    '''
    d = ((x2-x1)**2 + (y2-y1)**2)**0.5
    return d

