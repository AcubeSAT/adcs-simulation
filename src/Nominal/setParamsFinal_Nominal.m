%% ======= Satellite ========
dt = .1; %Timestep for Orbit Propagator
orbits=1;
orbitPeriod=5545;
% tf = orbits*orbitPeriod+0.9; %Total Simulation Seconds
tf =(5545+0.9); %Total Simulation Seconds
np = uint32((tf+dt)/dt); %Number of timesteps
q_desired = [ 1 0 0 0 ] ; %Desired quaternion
Q0 = [0.5; 0.5; 0.5; 0.5]; %Initial Quaternion in ECI frame
Q0 = Q0/norm(Q0); %Normalised Quaternion
init_bias = [.01;0.15;-.08]; % bias initialization
% vRot0 = [pi/4; pi/2; pi/8];
vRot0 = [0.035; 0.035; 0.035]; %Initial Angular Velocities from Body to ECI frame expressed in Body.
x0 = [Q0;vRot0]; %Initial state consists of [Quaternion;Angular Velocity]
% Random Initial State Estimation
Q0_hat = [.6;.1;-.7;.01];
Q0_hat = Q0_hat / norm(Q0_hat);
%vRot0_hat = [.7;-.8;.9];
bias_hat = [0;0;0];
x0_hat = [Q0_hat;bias_hat];

kd=1e-5; %Atrificial Detumbling gain

% torq = randn([3 np]).*1e-5; %Torque
% torq = randn([3 np]).*0; %Random Initial Torque
torq = zeros(3,1); %Random Initial Torque
magn_ref = 1.0e+04 * [2; 0.15; 3.25]; %Random Mag_Field Initial value
sun_ref = 1.0e1 * [.8;.2;.45];  %Random Sun_pos Initial value
model = SatelliteModel(dt, Const.I); %Initialize Satellite Model Class
real_model = real_SatelliteModel(dt, Const.I);
disturbancesEnabled = "total";
%% ======= Simulation Constants ========
n_dim = length(x0); % state length
n_dim_error = 6; % error state length
n_msr = 9; %Measurements length
Kp=100;
mtq_max = 0.2;
setDisturbances = "total";   % Set which disturbances you want to activate: tau_g, tau_ad, tau_sp, tau_rm, total, zero
rng(1); % Fix the random number generator for reproducible results
plotter_step=10;

%% ======= Orbit Propagation ========
satrec = orbit_init();
[xsat_ecf, vsat_ecf,xsat_eci,vsat_eci, sat_llh,eclipse, mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit, sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,satrec,argpm,nodem,inclm,mm,xnode,xinc] = orbit_sgp4(satrec,dt,tf+dt);

for(i=1:length(xsat_eci))
    xsat_eci_normalized(:,i) = xsat_eci(:,i)/norm(xsat_eci(:,i));
end

%eclipse = zeros(1,55460);
% save('mag_orbit_10.mat','mag_field_orbit');

%% ======= Albedo ========
albedo = load("sso_albedo.mat");  % Choose albedo depending on orbit
%albedo = load("iss_albedo.mat");

albedo = albedo.new;

%% ======= Kalman filter params ========
% Variances
Q = 0.5e-05*eye(n_dim_error,n_dim_error); % Variance of the process noise w[k]

%MGN noise 1e-3 (norm) | GYRO noise 1.57e-2| SUN noise 8.7e-3(norm)
% R_coeff=[1e-6;1e-6;1e-6;5e-5;5e-5;5e-5;1.2e-5;1.2e-5;1.2e-5];   
R_coeff=[1.83e-6;1.83e-6;1.83e-6;0;0;0;0;0;0];   
R = R_coeff.*eye(n_msr,n_msr); % Variance of the measurement noise v[k]

% Gyro bias std dev
sigma_u = 7.7570e-04;
% Gyro white noise std dev
sigma_v = 0.0026;

% R Variances used in EKF
% R_hat_coeff=[1e-3;1e-3;1e-3;8e-3;8e-3;8e-3;5e-3;5e-3;5e-3];
R_hat_coeff=[.5e-3;.5e-3;.5e-3;4e-3;4e-3;4e-3;1e-3;1e-3;1e-3];
R_hat = R_hat_coeff.*eye(n_msr,n_msr);

% Initialize Covariance matrix
P0 = 1*eye(n_dim_error,n_dim_error);

use_analytic_jacob = true;



