## Orbit Propagator
This folder includes the SGP4 Orbit Propagator software, directly cloned from the [Aerodynamic Software by David Vallado](https://celestrak.com/software/vallado-sw.php). This model is used to estimate the predicted satellite position on the Earth-Centered Inertial frame based on the satellite time and TLE file.

#### The used TLE files can be found in SGP4/DemoData

Important functions:
#### SGP4/twoline2rvMOD.m 
Which converts the two line element set character string data to variables and initializes the sgp4 variables.

#### SPG4/sgp4_2.m 
Which is a modified alternative to the sgp4.m function. This procedure is the sgp4 prediction model from space command. The mentioned modification allows the return of orbital parameters, as well as satellite position and velocity vectors.
