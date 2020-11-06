%% ======= Satellite ========
dt = 1; %Timestep for Controller & Orbit Propagator
dt_model= .005; %Timestep for Model
orbits=8;
orbitPeriod=5545;
tf = orbits*orbitPeriod+0.9; %Total Simulation Seconds
np = uint32((tf+dt)/dt); %Number of timesteps
Q0 = [0.3757;0.5983;-0.2364;-0.6671]; %Initial Quaternion in ECI frame
Q0 = Q0/norm(Q0); %Normalised Quaternion
vRot0 = [pi/2; pi/2; pi/8];
%vRot0 = [pi/6; -pi/6; pi/6]; %Initial Angular Velocities from Body to ECI frame expressed in Body.
x0 = [Q0; vRot0];

% torq_file = load('TorqueDet_6000_01.mat','T_magnetic');
% torq=torq_file.T_magnetic;
torq = randn([3 np]).*0; %Random Initial Torque
magn_ref = 1.0e+04 * [2; 0.15; 3.25]; %Random Mag_Field Initial value
sun_ref = 1.0e1 * [.8;.2;.45];  %Random Sun_pos Initial value
% model = SatelliteModel(dt, Const.I); %Initialize Satellite Model Class
real_model = real_SatelliteModel_Bdot(dt_model, Const.I);

%% ======= Simulation Constants ========
n_msr = 3; %Measurements length

setDisturbances = "total";   % Set which disturbances you want to activate: tau_g, tau_ad, tau_sp, tau_rm, total, zero
rng(1); % Fix the random number generator for reproducible results
plotter_step=10;

%% ======= Orbit Propagation ========
satrec = orbit_init();
% [xsat_ecf, vsat_ecf,xsat_eci,vsat_eci, sat_llh,eclipse, mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit, sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,satrec,argpm,nodem,inclm,mm,xnode,xinc] = orbit_sgp4(satrec,dt_model,tf+dt/dt_model);
[xsat_ecf, vsat_ecf,xsat_eci,vsat_eci, sat_llh,eclipse, mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit, sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,satrec,argpm,nodem,inclm,mm,xnode,xinc] = orbit_sgp4(satrec,dt,tf+dt);

%Kp = 10^5 * 4*pi / orbitPeriod * (1 + sin(inclm(1))) * Const.I(3,3);
Kp = 1;
% save('mag_orbit_10.mat','mag_field_orbit');

%% ======= Sensors ========
%MGN noise 1e-3 (norm) | GYRO noise 1.57e-2| SUN noise 8.7e-3(norm)
%MGN noise 3e-8
% R_coeff=[1e-6;1e-6;1e-6;5e-5;5e-5;5e-5;1.2e-5;1.2e-5;1.2e-5];   
R_coeff=[9e-16;9e-16;9e-16];
R = R_coeff.*eye(n_msr,n_msr); % Variance of the measurement noise v[k]

% Gyro bias std dev
sigma_u = 7.7570e-04;
% Gyro white noise std dev
sigma_v = 0.0026;