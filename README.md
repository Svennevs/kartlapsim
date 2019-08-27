# Go-kart lap time simulation
This repo contains a simple lap time simulation algorithm, applied to go-karts.

##Method
It's a quasi-steady-state method, vehicle model is reduced to a point-mass. Lateral acceleration limit is determined by tyre friction.
Longitudinal acceleration limit is determined by weight distribution, brake config (4 wheels or rear wheels) and tyre friction.
Combined acceleration limit determined by an ellipse, defined by long and lat acc limits. Other model features:
* Apex definition: intersection of ax for max ay with gradient of v for max ay.
* Forward/backward calculation from apexes: longitudinal acceleration at given speed and curvature includes slip angle drag, rolling and aerodynamic resistance and track inclination.
* Euler integration used in forward-backwards calculation
* Longitudinal tyre slip in presence of speed dependent engine torque limit (defined by torque curve) is solved using the brentq solver of scipy.optimize.

##Applications
* Example_1 is a simple lap time simulation for a mini kart at the circuit of Mariembourg
* Example_2 contains a scan of final drive ratio options for this circuit
The curvature data was obtained using gps logging (Mariembourg, 2019). The engine torque characteristic comes from dyno testing.

##Dependencies
Tested for python 3.6 (should work above)
---Packages:---
* matplotlib
* numpy
* math
* pickle (or cPickle for Python 2)
* scipy (specifically for scipy.optimize)
* itertools

##Develop branch
To improve precision of the method, it was attempted to include a two-track model.
Solving the apex solution and acceleration is now done numerically. So far the method was not succesfull.
With the scipy optimizers no reliable optima could be found.


