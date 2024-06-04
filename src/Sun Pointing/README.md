## Sun Pointing

The Sun Pointing simulation consists of three files:
- `main_Sun_Pointing`
- `setParamsFinal_Sun_Pointing`
- `PD_Sun_Pointing`

The main executable file is `main_Sun_Pointing`. In this file, there is the option to choose which diagrams will be plotted.

The important parameters of the Sun Pointing simulation are set in `setParamsFinal_Sun_Pointing`. Some of these parameters are:
- `orbits`: Number of orbits that the simulation will run
- `Q0`: The initial quaternion
- `vRot0`: The initial angular velocity
- `DisturbancesEnabled`: Which disturbances will be activated
- `init_bias`: Bias initialization
- `base_albedo`: Albedo depending on orbit 
- `sigma_u` & `sigma_v`: Bias standard deviation & White noise standard deviation 
- `R`, `R_hat` & `Q`: Variance matrices for Kalman filter

Parameters and constants that concern the satellite in general are set in `src/constants/constants`.

The file `PD_Sun_Pointing` is a function that implements the PD controller used to calculate the torque needed to achieve Sun Pointing. The gains `Kp_gain` and `Kd_gain` of the PD controller are set in this function.

The sun pointing simulation has the following dependencies:
- Actuator functions in `src/utils/system_model/actuators`
- Disturbances functions in `src/utils/system_model/disturbances`
- System model classes in `src/utils/system_model`
- Orbit functions in `src/orbit functions` and the orbit propagator in `src/orbit propagator`
- Kalman filter functions in `src/utils/kf_lib`
- Albedo toolbox in `src/Albedo Toolbox-1.0`
- Cookie functions used by the Kalman filter in `src/cookies`
- Gyro modeling function in `src/gyro modeling`
- Sun vector functions in `src/sun vectors`
- Function for solving Wahba's problem in `src/wahba`
- Mathematical functions in `src/utils/quaternions`, `src/utils/rotations`, and `src/transformations`

