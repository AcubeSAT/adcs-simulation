## Nadir Pointing
The Nadir Pointing simulation consists of three files:
- `main_Nadir_Pointing`
- `Nadir_Pointing_function`
- `setParamsFinal_Nadir_Pointing`
- `PD_Nadir_Pointing`
- `trigger_N2D`
- `Q_selection`
There is also the  subdirectory `Genetic Algorithm` that contains the files `fitness_function_Nadir_Pointing` and `Genetic_algorithm`. 

The main executable file is `main_nadir_Pointing`. The gains `Kp_gain` and `Kd_gain` of the PD controller are set in this file.

The function `Nadir_Pointing_function` implements the simulation loop. In this file, there is the option to choose which diagrams will be plotted.

The function `PD_Nadir_Pointing` implements the PD controller used to calculate the torque needed to achieve Nadir Pointing.

The important parameters of the Nadir Pointing simulation are set in `setParamsFinal_Nadir_Pointing`. Some of these parameters are:
- `orbits`: Number of orbits that the simulation will run
- `Q0`: The initial quaternion
- `vRot0`: The initial angular velocity
- `DisturbancesEnabled`: Which disturbances will be activated
- `init_bias`: Bias initialization
- `base_albedo`: Albedo depending on orbit 
- `sigma_u` & `sigma_v`: Bias standard deviation & White noise standard deviation 
- `R`, `R_hat` & `Q`: Variance matrices for Kalman filter

The function `trigger_N2D` calculates if Detumbling is ready to be activated. The angular velocity threshold for changing modes `N2D_threshold` is set in this file. 

The function `Q_selection` sets the Q and R variance matrices used by the MEKF.

#### Genetic Algorithm
This folder contains the files `fitness_function_Nadir_Pointing` and `Genetic_algorithm`. 

The `Genetic_algorithm` is used for the gain tuning of the PD controller for Nadir Pointing. 

The `fitness_function_Nadir_Pointing` is used by the `Genetic_algorithm`. It returns the fitness parameter, which is used to calculate the controller's gains.

