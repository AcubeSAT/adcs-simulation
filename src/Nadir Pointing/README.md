## Nadir Pointing
This folder includes all files simulating Nadir Pointing mode.

#### main_Nadir_Pointing
Executable file for Nadir Pointing simulation.

#### Nadir_Pointing_function
Main function implementing the simulation.

#### setParamsFinal_Nadir_Pointing
Script that defines all nadir pointing related parameters needed for the simulation.

#### PD_Nadir_Pointing
Function that implements the PD controller used by nadir pointing mode, which calculates the desired torque to be applied.

#### trigger_N2D
Funciton that calculates if Detumbling mode is ready to be activated, when in Nadir Pointing mode.

#### Q_selection
Function that sets Q and R matrices, uused by MEKF.

#### Genetic Algorithm
Folder that contains the files implementing the genetic algorithm used for tuning the PD controller's gains.