% ========================================================================
%   Script for parameters definition in Detimbling Mode
% ========================================================================

%% ======= Satellite ========
cycle_duration = 1;                                         % Timestep for Controller & Orbit Propagator
dt_model = .001;                                            % Timestep for Model
orbits = 3;
orbitPeriod = 5545;
Total_simulation_time = orbits * orbitPeriod;               % Total Simulation Seconds
Q0 = [0.3757; 0.5983; -0.2364; -0.6671];                    % Initial Quaternion in ECI frame
Q0 = Q0 / norm(Q0);                                         % Normalised Quaternion
vRot0 = [pi / 6; -pi / 6; pi / 6];                          % Initial Angular Velocities from Body to ECI frame expressed in Body.
x0 = [Q0; vRot0];
real_model = real_SatelliteModel_Bdot(dt_model, Const.I);   % Initialize Satellite Model Class

%% ======= Simulation Constants ========
setDisturbances = "total";  % Set which disturbances you want to activate: tau_g, tau_ad, tau_sp, tau_rm, total, zero
rng(1);                     % Fix the random number generator for reproducible results
plotter_step = 10;
%Kp = 10^5 * 4*pi / orbitPeriod * (1 + sin(inclm(1))) * Const.I(3,3);
Kp = 1;

%% ======= Orbit Propagation ========
satrec = orbit_init();
[xsat_ecf, vsat_ecf, xsat_eci, vsat_eci, sat_llh, eclipse, mag_field_ned, mag_field_eci, mag_field_ecef, mag_field_orbit, sun_pos_ned, sun_pos_eci, sun_pos_ecef, sun_pos_orbit, satrec, argpm, nodem, inclm, mm, xnode, xinc] = orbit_sgp4(satrec, cycle_duration, cycle_duration+Total_simulation_time);

%% ======= Sensors ========
n_msr = 3; % Measurements length
R_coeff = [9e-16; 9e-16; 9e-16];
R = R_coeff .* eye(n_msr, n_msr); % Variance of the measurement noise v[k]


%% Parameters for threshold limits

  total_limit = 520;
  exceptions_limit = 30;
