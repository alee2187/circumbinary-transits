# circumbinary-transits

## Overview

This repository contains some relatively simple code written in Python to simulate circumbinary orbits and track transits. The objective of this code is to take a circumbinary system and simulate the patterns in the transits to look for any irregularities within the transits. By simulating these orbits, we hope to be able to create some filters or signatures that can be used to sift through existing telescope data to hunt for transiting circumbinary planets. 

The results of this search have not been thoroughly investigated yet, but once those are available, they will be available on this GitHub page. 

## Scripts

### functions.py

This file contains the functions that are used in the simulation. 

### system.py

This file houses the code for the System class, which is used to model a circumbinary planet system. It is instantiated by specifying the mass ratio, binary and orbital eccentricity, radii of each body, and the inclination. The simulation runs for a default of three periods, and takes 10^5 data points per simulation. There is also functionality to plot the orbits of the system in 3D. 