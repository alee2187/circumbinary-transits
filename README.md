# circumbinary-transits

## Overview

This repository contains some relatively simple code written in Python to simulate circumbinary orbits and track transits. The objective of this code is to take a circumbinary system and simulate the lightcurve of the system to potentially detect some patterns in the transits. By simulating these orbits, we hope to be able to create some filters or signatures that can be used to sift through existing telescope data to hunt for transiting circumbinary planets. 

The results of this search have not been thoroughly investigated yet, but once those are available, they will be available on this GitHub page. 

## Scripts

### constants.py

This file contains some constants used in the other files (e.g., fundamental constants, conversions between units). 

### functions.py

This file contains the functions that are used in the main code (e.g., time() generates an array of true anomaly as a function of time, HW99() computes the Holam-Wiegert 1999 critical semimajor axis). 

### simulation.py

This file houses the main simulation, which takes a given circumbinary system and prints the number of transits of the planet over each star per epoch. The code is set to three epochs, but this can be changed. There is also functionality to graph the orbits of the system in 3D. 

### kepler16b.txt

This file contains the data for the Kepler-16b circumbinary planet, which gets imported into the simulation. 

## Example

The simulation.py file contains an example of the signatures being produced for the Kepler-16b system. Running the code gives us that the planet transits star A once per epoch, but transits star B 3, 0, and 1 time(s) in each of the three epochs, respectively. 
