# circumbinary-transits

## Overview

This repository contains some relatively simple code written in Python to simulate circumbinary orbits and track transits. The objective of this code is to take a circumbinary system and simulate the patterns in the transits to look for any irregularities within the transits. By simulating these orbits, we hope to be able to create some filters or signatures that can be used to sift through existing telescope data to hunt for transiting circumbinary planets. 

The results of this search have not been thoroughly investigated yet, but once those are available, they will be available on this GitHub page. 

## Scripts

### functions.py

This file contains the functions that are used in the simulation. 

### system.py

This file houses the code for the <code>System</code> class, which is used to model a circumbinary planet system. It is instantiated by specifying the mass ratio, binary and orbital eccentricity, radii of each body, and the inclination. The simulation runs for a default of three periods, and takes 10<sup>5</sup> data points per simulation. There is also functionality to plot the orbits of the system in 3D. 

### test.py

This file contains an example system, with numbers matching the circumbinary system Kepler 16b. 
<!-- CHANGE THE NUMBERS IN TEST.PY TO MATCH KEPLER -->

## Example

Here are some example plots of the Kepler 16b system, the outputs of <code>test.py</code>. 
<!-- work in progress. need to plot and show some nice graphics. maybe write a gif to show the orbit in real time around a transit? -->
