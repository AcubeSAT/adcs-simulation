% ========================================================================
%   Script for parameters definition in Z-Thomson Mode
% ========================================================================

function Param = setParams_ZThomson(I)

%% ======= Satellite ========
cycle_duration = 0.1;                                         % Timestep for Controller & Orbit Propagator
dt = 0.1;
dt_model = .001;                                            % Timestep for Model
orbits = 1;
orbitPeriod = 5545;
tf = orbits * orbitPeriod;  %Total Simulation Seconds
Total_simulation_time = orbits * orbitPeriod;               % Total Simulation Seconds
Q0 = [0.3757; 0.5983; -0.2364; -0.6671];                    % Initial Quaternion in ECI frame
Q0 = Q0 / norm(Q0);                                         % Normalised Quaternion
vRot0 = [0.05 ; 0.05 ; 0.016];                          % Initial Angular Velocities from Body to ECI frame expressed in Body.
x0 = [Q0; vRot0];
%real_model = real_SatelliteModel_Bdot(dt_model, I);   % Initialize Satellite Model Class                          % Desired angular velocity [rad/s]
w_ref = 0.087;

%% ======= Orbit Propagation ========
satrec = orbit_init();
[xsat_ecf, vsat_ecf, xsat_eci, vsat_eci, sat_llh, eclipse, mag_field_ned, mag_field_eci, mag_field_ecef, mag_field_orbit, sun_pos_ned, sun_pos_eci, sun_pos_ecef, sun_pos_orbit, satrec, argpm, nodem, inclm, mm, xnode, xinc] = orbit_sgp4(satrec,dt,tf+dt);

%% ============== Calculations ====================================
init_bias = [.01; 0.15; -.08];      % bias initialization
vRot0 = [0; 0; 0];                  %Initial Angular Velocities from Body to ECI frame expressed in Body.
x0 = [Q0; vRot0];                   %Initial state consists of [Quaternion;Angular Velocity]
Q0_hat = [.6; .1; -.7; .01];        %Random Initial State Estimation
Q0_hat = Q0_hat / norm(Q0_hat);
bias_hat = [0; 0; 0];
x0_hat = [Q0_hat; bias_hat];
model = SatelliteModel(dt, I);      %Initialize Satellite Model Class
real_model = real_SatelliteModel_Bdot(dt_model, I);
N_Timesteps = 10;

%% ======= Simulation Constants ========
setDisturbances = "total";  % Set which disturbances you want to activate: tau_g, tau_ad, tau_sp, tau_rm, total, zero
rng(1);                     % Fix the random number generator for reproducible results
plotter_step = 10;
%Kd = 10^5 * 4*pi / orbitPeriod * (1 + sin(inclm(1))) * Const.I(3,3);
Kd = 100;
Ks = 1;

sigma_u = 3.4434e-4;
sigma_v = 0.0004;

%% ======= Albedo ========
%albedo = load("sso_albedo.mat");  % Choose albedo depending on orbit
%albedo = load("iss_albedo.mat");
base_albedo = load("SSO_500_6PM_1_Orbit_051231.mat");
base_albedo = base_albedo.new;
base_albedo_inaccurate = load("SSO_500_6PM_1_Orbit_050101");
base_albedo_inaccurate = base_albedo_inaccurate.new;
albedo = base_albedo;
albedo_inaccurate = base_albedo_inaccurate;
if orbits > 1
    for i = 1:orbits - 1
        albedo = [albedo, base_albedo];
        albedo_inaccurate = [albedo_inaccurate, base_albedo_inaccurate];
    end
end
%albedo_inaccurate = load("sso_albedo.mat");


%% ======= Sensors ========
n_msr = 3; % Measurements length
R_coeff_old = [9e-16; 9e-16; 9e-16];
R_old = R_coeff_old .* eye(n_msr, n_msr); % Variance of the measurement noise v[k]

%% ======= Kalman filter params ========
Q = 1e-4 * diag([1, 1, 1, 1e-3, 1e-3, 1e-3]);

%MGN noise 1e-3 (norm) | GYRO noise 1.57e-2| SUN noise 8.7e-3(norm)
% R_coeff=[1e-6;1e-6;1e-6;5e-5;5e-5;5e-5;1.2e-5;1.2e-5;1.2e-5];
R_coeff = [1.83e-6; 1.83e-6; 1.83e-6; 0; 0; 0; 0; 0; 0];
number_of_measurements = 9; %Measurements length
R = R_coeff .* eye(number_of_measurements, number_of_measurements); % Variance of the measurement noise v[k]

% R Variances used in MEKF
% R_hat_coeff=[1e-3;1e-3;1e-3;8e-3;8e-3;8e-3;5e-3;5e-3;5e-3];
R_hat_coeff = [.5e-3; .5e-3; .5e-3; 1e-3; 1e-3; 1e-3];
R_hat = R_hat_coeff .* eye(6, 6);
% Initialize Covariance matrix
n_dim_error = 6; % error state length
P0 = 1 * eye(n_dim_error, n_dim_error);
use_analytic_jacob = true;

%% Parameters for threshold limits

total_limit = 520;
exceptions_limit = 30;

%%  Passing the values of the parameters in a struct.
Param.N_Timesteps = N_Timesteps;
Param.tf = tf;
Param.w_ref = w_ref;
Param.Kd = Kd;
Param.Ks = Ks;
Param.satrec = satrec;
Param.R_old = R_old;
Param.tst = Total_simulation_time;
Param.cycle_duration = cycle_duration;
Param.dt = dt;
Param.orbits = orbits;
Param.x0 = x0;
Param.x0_hat = x0_hat;
Param.real_model = real_model;
Param.model = model;
Param.albedo = albedo;
Param.albedo_inaccurate = albedo_inaccurate;
Param.setDisturbances = setDisturbances;
Param.eclipse = eclipse;
Param.mag_field_eci = mag_field_eci;
Param.mag_field_orbit = mag_field_orbit;
Param.sun_pos_eci = sun_pos_eci;
Param.sun_pos_orbit = sun_pos_orbit;
Param.argpm = argpm;
Param.nodem = nodem;
Param.inclm = inclm;
Param.mm = mm;
Param.init_bias = init_bias;
Param.Q = Q;
Param.R_coeff = R_coeff;
Param.R = R;
Param.R_hat = R_hat;
Param.sigma_u = sigma_u;
Param.sigma_v = sigma_v;
Param.P0 = P0;
Param.number_of_measurements = number_of_measurements;
Param.use_analytic_jacob = use_analytic_jacob;
Param.xsat_eci = xsat_eci;
Param.total_limit= total_limit;


end
