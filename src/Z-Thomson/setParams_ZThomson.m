% ========================================================================
%   Function for parameters definition in Z Thomson simulation.
%
%   Inputs
%     I         - Inertia matrix
%   Outputs:
%     Param     - Struct containing simulation parameters
% ========================================================================

function Param = setParams_ZThomson(I)

    %% ======= Satellite ========

    dt = .1; %Timestep for Orbit Propagator
    orbits = 3;
    orbitPeriod = 5545;
    tf = orbits * orbitPeriod;  %Total Simulation Seconds
    q_desired = [1, 0, 0, 0];   %Desired quaternion
    N_Timesteps = 10;           % Number of timesteps per cycle

    Kd = 100;
    Ks = 1;
    w_ref = 0.087;

    %% ======= Orbit Propagation ========

    [satrec, ~] = orbit_init();
    [~, ~, xsat_eci, ~, ~, eclipse, ~, mag_field_eci, ~, mag_field_orbit, ~, sun_pos_eci, ~, sun_pos_orbit, ~, argpm, nodem, inclm, mm, ~, ~] = orbit_sgp4(satrec, dt, tf+dt);
    %[~,~,xsat_eci,~,~,eclipse,~,mag_field_eci,~,mag_field_orbit,~,sun_pos_eci,~,sun_pos_orbit,~,argpm,nodem,inclm,mm,~,~] = orbit_sgp4_offset(satrec,dt,tf+dt,1000);




    %% ======================= Testing Initial Quaternions =============================
    %Q0 = [0.5; -0.5; 0.5; 0.5];
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

    %% =============== Desired Quaternion =================================

    Q0 = [-0.4493; 0.1189; -0.8854; 0.0122]; % [1, 0, 0, 0] TLE -> 6PM 2-offset %Initial Quaternion in ECI frame

    %% ============== Calculations ====================================

    Q0 = Q0 / norm(Q0);                     %Normalised Quaternion
    init_bias = [.01; 0.15; -.08];          % bias initialization
    vRot0 = [0.05; 0.05; 0.05];                      %Initial Angular Velocities from Body to ECI frame expressed in Body.
    x0 = [Q0; vRot0];                       %Initial state consists of [Quaternion;Angular Velocity]
    Q0_hat = [.6; .1; -.7; .01];            %Random Initial State Estimation
    Q0_hat = Q0_hat / norm(Q0_hat);
    bias_hat = [0; 0; 0];
    x0_hat = [Q0_hat; bias_hat];
    model = SatelliteModel(dt, I);          %Initialize Satellite Model Class
    real_model = real_SatelliteModel(dt, I);

    %% ======= Simulation Constants ========

    disturbancesEnabled = "total";      % Set which disturbances you want to activate: tau_g, tau_ad, tau_sp, tau_rm, total, zero
    rng(1);                             % Fix the random number generator for reproducible results
    
    %Gyro Noise Parameters (SCHA63T)
    ARW=2.27e-5; %(in (rad/sec)/sqrt(Hz))
    RRW=3e-9;  %((in (rad/sec)*sqrt(Hz)))
    BI=2.8e-7; %(in rad/sec)

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

    %%  Parameters for threshold limits
    total_limit = 20;
    exceptions_limit = 2;

    %%  Passing the values of the parameters in a struct.
    Param.satrec=satrec;
    Param.dt = dt;
    Param.orbits = orbits;
    Param.tf = tf;
    Param.q_desired = q_desired;
    Param.x0 = x0;
    Param.x0_hat = x0_hat;
    Param.real_model = real_model;
    Param.model = model;
    Param.albedo = albedo;
    Param.albedo_inaccurate = albedo_inaccurate;
    Param.disturbancesEnabled = disturbancesEnabled;
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
    Param.P0 = P0;
    Param.number_of_measurements = number_of_measurements;
    Param.use_analytic_jacob = use_analytic_jacob;
    Param.xsat_eci = xsat_eci;
    Param.total_limit= total_limit;
    Param.exceptions_limit= exceptions_limit;
    Param.N_Timesteps= N_Timesteps;
    Param.Kd = Kd;
    Param.Ks = Ks;
    Param.w_ref = w_ref;
    Param.ARW = ARW;
    Param.RRW = RRW;
    Param.BI= BI;
end