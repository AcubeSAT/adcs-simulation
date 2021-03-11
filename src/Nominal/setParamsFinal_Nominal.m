function Param = setParamsFinal_Nominal(I)


%% ======= Satellite ========
dt = .1; %Timestep for Orbit Propagator
orbits=1;
orbitPeriod=5545;
tf = orbits*orbitPeriod; %Total Simulation Seconds
%tf =(5545+0.9); %Total Simulation Seconds
np = uint32((tf+dt)/dt); %Number of timesteps
q_desired = [ 1 0 0 0 ] ; %Desired quaternion

%% ======= Orbit Propagation ========
[satrec, x] = orbit_init();
[xsat_ecf, vsat_ecf,xsat_eci,vsat_eci, sat_llh,eclipse, mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit, sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,satrec,argpm,nodem,inclm,mm,xnode,xinc] = orbit_sgp4(satrec,dt,tf+dt);

for(i=1:length(xsat_eci))
    xsat_eci_normalized(:,i) = xsat_eci(:,i)/norm(xsat_eci(:,i));
end

%eclipse = zeros(1,55460);
% save('mag_orbit_10.mat','mag_field_orbit');

%% ======================= Testing Initial Quaternions =============================
% if x >= 0 && x <= 0.25
%Q0 = [-0.4493; 0.1189; -0.8854; 0.0122];% [1, 0, 0, 0] TLE -> 6PM 2-offset
% elseif x > 0.25 && x <= 0.5
%     Q0 = [-0.4371; 0.3440; -0.8244; -0.1045];% [1, 0, 0, 0] TLE -> 8PM 2-offset
% elseif x > 0.5 && x <= 0.75
%     Q0 = [-0.4197; 0.4487; -0.7725; -0.1607];% [1, 0, 0, 0] TLE -> 9PM 2-offset  
% else
%     Q0 = [0.3638; -0.6333; 0.6300; 0.2639];% [1, 0, 0, 0] TLE -> 11PM 2-offset
% end
%Q0 = [0.5; -0.5; 0.5; 0.5]; %Initial Quaternion in ECI frame
%Q0 = [0.0894; -0.0058; 0.0641; -0.9939];% [1, 0, 0, 0] TLE -> 29_08
%Q0 = [0.9928; -0.0641; -0.0054; 0.1007];% [1, 0, 0, 0] TLE -> 6PM 0-offset
%Q0 = [-0.4493; 0.1189; -0.8854; 0.0122];% [1, 0, 0, 0] TLE -> 6PM 2-offset
%Q0 = [0.0619; 0.3476; -0.9310; 0.0920];% [1, 0, 0, 0] TLE -> 8PM 3-offset
%Q0 = [-0.4371; 0.3440; -0.8244; -0.1045];% [1, 0, 0, 0] TLE -> 8PM 2-offset
%Q0 = [0.0493; 0.4662; -0.8777; 0.0993];% [1, 0, 0, 0] TLE -> 9PM 3-offset
%Q0 = [-0.4197; 0.4487; -0.7725; -0.1607];% [1, 0, 0, 0] TLE -> 9PM 2-offset
%Q0 = [0.0219; 0.6775; -0.7271; 0.1087];% [1, 0, 0, 0] TLE -> 11PM 3-offset
%Q0 = [0.3638; -0.6333; 0.6300; 0.2639];% [1, 0, 0, 0] TLE -> 11PM 2-offset
%Q0 = [0.8315; -0.6103; -0.0033; 0.8585];

%% =============== Random Initial Quaternion ==================================
%Q0 = [0.05; 0.45; -0.88; 0.1]; %-> Adaptive good
%Q0 = [0.4; 0.3; -0.7; 0.5]; %-> Non-adaptive good

%% =============== Desired Quaternion =================================
Q0 = [-0.4493; 0.1189; -0.8854; 0.0122];% [1, 0, 0, 0] TLE -> 6PM 2-offset

%% ============== Calculations ====================================

Q0 = Q0/norm(Q0); %Normalised Quaternion
init_bias = [.01;0.15;-.08]; % bias initialization
%vRot0 = [0; 0; 0];
%vRot0 = [pi/4; pi/2; pi/8];
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
model = SatelliteModel(dt, I); %Initialize Satellite Model Class
real_model = real_SatelliteModel(dt, I);
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



%% ======= Albedo ========
%albedo = load("sso_albedo.mat");  % Choose albedo depending on orbit
%albedo = load("iss_albedo.mat");
base_albedo = load("SSO_500_6PM_1_Orbit_051231.mat");
base_albedo = base_albedo.new;
base_albedo_inaccurate = load("SSO_500_6PM_1_Orbit_050101");
base_albedo_inaccurate = base_albedo_inaccurate.new;
albedo = base_albedo;
albedo_inaccurate = base_albedo_inaccurate;
if orbits>1
    for i=1:orbits-1
        albedo = [albedo base_albedo];
        albedo_inaccurate = [albedo_inaccurate base_albedo_inaccurate];
    end
end


%albedo_inaccurate = load("sso_albedo.mat"); 


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

%%  Passing the values of the parameters in a struct.

    Param.dt = dt;
    Param.orbits = orbits;
    Param.orbitPeriod = orbitPeriod;
    Param.tf = tf;
    Param.np = np;
    Param.q_desired = q_desired;
    Param.x0 = x0;
    Param.x0_hat = x0_hat;
    Param.torq = torq;
    Param.magn_ref = magn_ref;
    Param.sun_ref = sun_ref;
    Param.real_model = real_model;
    Param.model = model;
    Param.n_msr = n_msr;
    Param.n_dim = n_dim;
    Param.kd = kd;
    Param.albedo = albedo;
    Param.albedo_inaccurate = albedo_inaccurate;
    Param.n_dim_error = n_dim_error;

    Param.disturbancesEnabled = disturbancesEnabled;
    Param.setDisturbances = setDisturbances;
    Param.plotter_step = plotter_step;
    Param.mtq_max = mtq_max;
    

    Param.xsat_ecf = xsat_ecf;
    Param.vsat_ecf = vsat_ecf;
    Param.xsat_eci = xsat_eci;
    Param.vsat_eci = vsat_eci;
    Param.sat_llh = sat_llh;
    Param.eclipse = eclipse;
    Param.mag_field_ned = mag_field_ned;
    Param.mag_field_eci = mag_field_eci;
    Param.mag_field_ecef = mag_field_ecef;
    Param.mag_field_orbit = mag_field_orbit;
    Param.sun_pos_ned = sun_pos_ned;
    Param.sun_pos_eci = sun_pos_eci;
    Param.sun_pos_ecef = sun_pos_ecef;
    Param.sun_pos_orbit = sun_pos_orbit;
    Param.satrec = satrec;
    Param.argpm = argpm;
    Param.nodem = nodem;
    Param.inclm = inclm;
    Param.mm = mm;
    Param.xnode = xnode;
    Param.xinc = xinc;

    Param.init_bias = init_bias;
    Param.Q = Q;
    Param.Kp = Kp;
    Param.R_coeff = R_coeff;
    Param.R = R;
    Param.R_hat = R_hat;
    Param.sigma_u = sigma_u;
    Param.sigma_v = sigma_v;
    Param.P0 = P0;
    Param.use_analytic_jacob = use_analytic_jacob;
    Param.sat_llh=sat_llh;


end



