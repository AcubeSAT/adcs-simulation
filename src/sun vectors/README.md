## Sun Functions (Sun Vectors)
This folder includes the Sun Position, Eclipse calculation and sun sensors modelling. The Sun Position and Eclipse models are used to estimate the predicted sun position and visibility on the Earth-Centered Inertial frame based on the satellite time and position. The sun sensor functions simulate the attitude and all corresponding errors of the coarse (and fine) sun sensors.
+ The Sun Position software is directly cloned from the [Aerodynamic Software by David Vallado](https://celestrak.com/software/vallado-sw.php).
+ The Eclipse Calculations software is developed based on [Fundamentals of Astrodynamics and Applications, 4th ed.](http://libgen.rs/book/index.php?md5=2EE909BE4133AF468EAD7E83B030B80D) and [Eclipse Prediction Algorithms for Low-Earth-Orbiting Satellites](https://ieeexplore.ieee.org/document/7967707).

#### calculate eclipse
This algorithm calculates the extistence of eclipse at any time during in-orbit. 

#### css compensation
This algorithm models the 6 coarse sun sensor units along with the albedo effect error induced in their measurements. Additionally this function includes the albedo compensation used in the MEKF in order to estimate tha satellite's state with greater accuracy. 

#### css noise
This algorithm models the 6 coarse sun sensor units along with the albedo effect error induced in their measurements. These measurements are used in the correct phase of the MEKF in order to correctly estimate the satellite's attitude.

#### sun fov plots
This algorithm models a fine sun sensor unit.

#### sun position
This algorithm is used to estimate the predicted sun position on the Earth-Centered Inertial frame based on the satellite time. This prediction is utilized in the MEKF correct step and is compared to the sensor measurements to estimate the satellite orientation.
