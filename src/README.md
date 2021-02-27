## Albedo Matrices
This folder includes Albedo effect matrices created by the corresponding defined TLEs. All matrices are produced by using the [Albedo Toolbox's](https://gitlab.com/acubesat/adcs/adcs-mekf/-/tree/GainTuning/src/AlbedoToolbox-1.0), [Default Albedo function](https://gitlab.com/acubesat/adcs/adcs-mekf/-/blob/GainTuning/src/AlbedoToolbox-1.0/albedo.m). The matrices are stored exclusively due to the high computational need of the Albedo Toolbox and in order to ease the simulation analysis process.

## Albedo Toolbox
This folder includes the Albedo Toolbox is directly cloned from the [Website of Dan D. V. Bhanderi](https://www.bhanderi.dk/downloads/). It is used to calculate the introduced albedo effect of the earth on the satellite body frame, which megatively affects the coarse sun sensor measurements.

## Bdot
This folder includes he AOCS simulation analysis main files regarding the AOCS Detumbling Mode.

The purpose of the AOCS Detumbling Mode is to autonomously dissipate and control the satellite's angular velocity. This is achieved with a B-dot controller directly using measurements from the high precision 3-axis magnetometer. In case of magnetometer failure, the redundant 3-axis magnetometer unit integrated in the iMTQ board is used. Actuation is performed solely by the magnetorquers. No satellite orbit propagation, initialized by a TLE file, is performed during this mode.

## Nominal
This folder includes the AOCS simulation analysis main files regarding the AOCS Nominal Mode while Nadir Pointing.

The purpose of the AOCS nominal mode is nadir-pointing according to the requirements. Nadir-pointing enables the reliable and continuous data transmission via the S-band to the GS. Nadir-pointing is achieved by determining and controlling the satellite's attitude. Attitude is determined by combining sensor measurements from the magnetometer, the sun sensors and the gyroscope. Attitude control is carried out by the magnetorquers and the reaction wheel.

## Sun Pointing (Nominal_Moustakas+Pappa)
This folder includes the AOCS simulation analysis main files regarding the AOCS Nominal Mode while Sun Pointing.

This is an optional attitude profile of AOCS nominal mode being utilized in specified orbits, for which the power budget does not meet the mission's requirements. The attitude determination and orbit propagator algorithms remain the same as in N-01(Nadir Pointing). A PD controller that aims to achieve sun-pointing is implemented and the commanded torque is split between both the actuators.

## Constants
This folder includes the AcubeSAT constants required in order to correctly simulate the satellite attitude during in-orbit conditions.

## Cookies
This folder includes functions that return structs in order to store all required input parameters for specific functions throughout the repository. The mentioned functions are mainly used to input parameters in the [Attitude Determination Filter (MEKF)](https://gitlab.com/acubesat/adcs/adcs-mekf/-/tree/GainTuning/src/utils/kf_lib).

## Gyro Modeling
This folder includes the gyroscope modeling functions, in order to correctly simulate its attitude and all corresponding errors.

## IGRF
This folder includes the International Geomagnetic Reference Field model, directly cloned from the [International Geomagnetic Reference Field (IGRF) Model by Drew Compston](https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field-igrf-model). This model is used to estimate the predicted magnetic field values on the Earth-Centered Inertial frame based on the satellite position.

## Orbit Functions
This folder includes the required functions in order to correctly interface the orbit propagator, IGRF, sun position and eclipse software whilst the satellite is in-orbit. The specified software are used to estimate the predicted reference models on the Earth-Centered Inertial frame based on the satellite time and TLE file.

## Orbit Propagator
This folder includes the SGP4 Orbit Propagator software, directly cloned from the [Aerodynamic Software by David Vallado](https://celestrak.com/software/vallado-sw.php). This model is used to estimate the predicted satellite position on the Earth-Centered Inertial frame based on the satellite time and TLE file.

## Sun Functions (Sun Vectors)
This folder includes the Sun Position, Eclipse calculation and sun sensors modelling. The Sun Position and Eclipse models are used to estimate the predicted sun position and visibility on the Earth-Centered Inertial frame based on the satellite time and position. The sun sensor functions simulate the attitude and all corresponding errors of the corse (and fine) sun sensors.
+ The Sun Position software is directly cloned from the [Aerodynamic Software by David Vallado](https://celestrak.com/software/vallado-sw.php).
+ The Eclipse Calculations software is developed based on [Fundamentals of Astrodynamics and Applications, 4th ed.](http://libgen.rs/book/index.php?md5=2EE909BE4133AF468EAD7E83B030B80D) and [Eclipse Prediction Algorithms for Low-Earth-Orbiting Satellites](https://ieeexplore.ieee.org/document/7967707).

## Transformations
This folder includes the required transformation rotation matrices between all defined frames with regards to the AcubeSAT AOCS design. The defined frames are:
+ Earth-Centered Inertial (ECI)
+ Earth-Centered, Earth-Fixed (ECEF/ECF)
+ local North-East-Down (NED)
+ Orbit frame

## Utils
This folder includes the utilities used for the AOCS simulation analysis, including Attitude Determination Kalman Filters, Kinematic & Dynamic Models and used Quaternion-based mathematics.

## Wahba
This folder includes the software solving Wahba's problem for Sensor Fusion used for Mathematican, and non-filtered, Attitude Determination. The software is used to estimate the initial satellite attitude when entering AOCS Nominal Mode. This estimastion is used to initialize the [Multiplicative Extended Kalman Filter (MEKF)](https://gitlab.com/acubesat/adcs/adcs-mekf/-/tree/GainTuning/src/utils/kf_lib/%40EKF), which is enabled and used subsequently.
