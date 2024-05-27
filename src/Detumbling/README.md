##  Detumbling

The detumbling simulation consists of three files:
- `main_Detumbling`
- `setParamsFinal_Detumbling`
- `trigger_D2N`

The main executable file is `main_Detumbling`. In this file, there is an option to choose which diagrams will be plotted.

The important parameters of the detumbling simulation are set in `setParamsFinal_Detumbling`. Some of these parameters are:
- `orbits`: Number of orbits that the simulation will run
- `Q0`: The initial quaternion
- `vRot0`: The initial angular velocity
- `setDisturbances`: Which disturbances will be activated

Parameters and constants that concern the satellite in general are set in `src/constants/constants`.

The file `trigger_D2N` is a function that calculates if Nadir Pointing is ready to be activated.

The detumbling simulation has the following dependencies:
- Actuator functions in `src/utils/system_model/actuators`
- Disturbances functions in `src/utils/system_model/disturbances`
- System model classes in `src/utils/system_model`
- Orbit functions in `src/orbit functions` and the orbit propagator in `src/orbit propagator`
- Mathematical functions in `src/utils/quaternions`, `src/utils/rotations`, and `src/transformations`