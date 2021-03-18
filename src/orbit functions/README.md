## Orbit Functions
This folder includes the required functions in order to correctly interface the orbit propagator, IGRF, sun position and eclipse software whilst the satellite is in-orbit. The specified software are used to estimate the predicted reference models on the Earth-Centered Inertial frame based on the satellite time and TLE file.

#### orbit_init
This function initializes the SGP4 orbit propagator based on the input TLE file.

#### orbit_main
This script is used to test the correct function of the SGP4 orbit propagator. It is not required for the general ADCS simulation analysis.

#### orbit_sgp4
This function propagates the satellite by utilizing the SPG4 orbit propagator software. Based on the propagator, the function returns magnetic field, sun position and eclipse calculations.

#### orbit_sgp4_offset
This function propagates the satellite by utilizing the SPG4 orbit propagator software. Based on the propagator, the function returns magnetic field, sun position and eclipse calculations. This function is based on the orbit_sgp4 function but can use offset time as an input variable

#### reference_vectors_calc
This function calculates the magnetic field and sun position in every ADCS defined frame, for every orbit propagator time and position pair.
