function [instant_error_perform, Time] = Nominal_Simulation(Kp_gain, Kd_gain) 

close all;
clc;

%clear variables;

%% Initialize Parameters Script
Const=constants();
Param = setParamsFinal_Nominal(Const.I);
%setParamsFinal_Nominal();
    dt = Param.dt;
    orbits = Param.orbits;
    %orbitPeriod = Param.orbitPeriod;
    tf = Param.tf;
    %np = Param.np;
    q_desired = Param.q_desired;
    x0 = Param.x0;
    x0_hat = Param.x0_hat;
    torq = Param.torq;
    %magn_ref = Param.magn_ref;
    %sun_ref = Param.sun_ref;
    real_model = Param.real_model;
    model = Param.model;
    n_msr = Param.n_msr;
    %n_dim = Param.n_dim;
    kd = Param.kd;
    albedo = Param.albedo;
    albedo_inaccurate = Param.albedo_inaccurate;

    disturbancesEnabled = Param.disturbancesEnabled;
    %setDisturbances = Param.setDisturbances;
    %plotter_step = Param.plotter_step;

    %xsat_ecf = Param.xsat_ecf;
    %vsat_ecf = Param.vsat_ecf;
    xsat_eci = Param.xsat_eci;
    %vsat_eci = Param.vsat_eci;
    %sat_llh = Param.sat_llh;
    eclipse = Param.eclipse;
    %mag_field_ned = Param.mag_field_ned;
    mag_field_eci = Param.mag_field_eci;
    %mag_field_ecef = Param.mag_field_ecef;
    mag_field_orbit = Param.mag_field_orbit;
    %sun_pos_ned = Param.sun_pos_ned;
    sun_pos_eci = Param.sun_pos_eci;
    %sun_pos_ecef = Param.sun_pos_ecef;
    sun_pos_orbit = Param.sun_pos_orbit;
    %satrec = Param.satrec;
    argpm = Param.argpm;
    nodem = Param.nodem;
    inclm = Param.inclm;
    mm = Param.mm;
    %xnode = Param.xnode;
    %xinc = Param.xinc;
    n_dim_error = Param.n_dim_error;
    sat_llh=Param.sat_llh;

    init_bias = Param.init_bias;
    Q = Param.Q;
    %Kp = Param.Kp;
    %R_coeff = Param.R_coeff;
    R = Param.R;
    R_hat = Param.R_hat;
    sigma_u = Param.sigma_u;
    sigma_v = Param.sigma_v;
    P0 = Param.P0;
    use_analytic_jacob = Param.use_analytic_jacob;
    
eclipse_data=0;

%% Initialize Global Parameters
global R_coils; R_coils=Const.R_coils;
global N_coils; N_coils=Const.N_coils;
global A_coils; A_coils=Const.A_coils;
global mt; mt = Const.mt;
global Jw; Jw=Const.Jw;
global Km; Km=Const.Km;
global Kv; Kv=Const.Kv;
global A; A=Const.A;
global b_friction; b_friction=Const.b_friction;
global c_friction; c_friction=Const.c_friction;
global Rb; Rb=Const.Rb;
global flag; flag = 0;
global Max_RW_torq; Max_RW_torq=Const.rw_max_torque;

%% Construct the EKF

Q_struct = load('idealQ.mat', 'IDEALQ');
Q_eclipse_load = Q_struct.IDEALQ;
n_params = length(x0_hat); %Init number of parameters
ekf = EKF(n_params, n_msr, @model.stateTransFun, @model.msrFun); %Init EKF
ekf.theta = x0_hat; %Init state estimation
ekf.P = P0; %Init Covariance matrix
ekf.setProcessNoiseCov(Q); %Q variance matrix
ekf.setMeasureNoiseCov(R_hat); %R variance matrix
ekf.setFadingMemoryCoeff(1.00); %This parameter defined the "memory" of the filter, 1 meaning regular
ekf.setPartDerivStep(0.001); %arithmetical jacobian step
if (use_analytic_jacob) %if analytical jacobian is used, define the functions
    ekf.setStateTransFunJacob(@model.stateTransFunJacob);
    ekf.setMsrFunJacob(@model.msrFunJacob);
end



%% Initialize simulation parameters

Time = 0:dt:tf;
x_real = zeros(7,length(Time)); % Real state
q_ob_data = zeros(4,length(Time));
x_hat_data = []; % Estimated state
x = x0(1:7);
x_real(:,1) = x0(1:7);
bias_data = init_bias;
gyro_noise_data = zeros(3,1);
t = 0;
n_steps = floor(length(Time)/10); % Number of time steps
timeflag_dz = 0;
init_AngVel_dz = 0;
init_accel_dz = 0;
rw_ang_momentum=0;
rw_ang_vel_rpm = zeros(1, length(Time));    % RW angular velocity in rpm
rw_accel = zeros(1, length(Time));          % RW acceleration
tau_mtq = zeros(3, length(Time));           % Torques produced by MTQs
tau_rw = zeros(1, length(Time));            % Torques produced by the RW
M_data = zeros(3, length(Time));
tau_dist = zeros(3, length(Time));
lambda=1;
tau_ad = zeros(3, length(Time));
tau_rm = zeros(3, length(Time));
tau_g = zeros(3, length(Time));
tau_sp = zeros(3, length(Time));
plotter_step = 1;
reps = length(Time);

%% Initialize matrices
AngVel_rw_radps = zeros(3, 1); %1 = old, 2 = cur, 3 = next
AngVel_rw_rpm = zeros(3, 1);
acceleration_rw = zeros(3, 1);
rm = [];

%% Next we initialize the bias estimation by solving Wahba's problem n times. 

bias_init_counter = 0;  % how many Wahba's have we solved for bias init
bias_wahba_loops = 2;   % Total times to be solved
quat_pos = zeros(4,bias_wahba_loops); % Wahba results are stored here
real_bias=init_bias;

entered_eclipse = false;
exited_eclipse = false;


for l=1:n_steps
    
    % Initial state when we get sun sensor measurement
    if eclipse((l-1)*10+1)==0
        
        %% Measurements
        y_real = real_model.msrFun(x_real(:,(l-1)*10+1),msrCookieFinal(mag_field_eci(:,(l-1)*10+1),...
            sun_pos_eci(:,(l-1)*10+1),eclipse((l-1)*10+1),[0;0;0]));
        y_noise = y_real + sqrt(R)*randn(size(y_real));
        [gyro_noise,real_bias] = gyro_noise_func(real_bias,dt,sigma_u,sigma_v);
        bias_data = [bias_data real_bias];
        gyro_noise_data = [gyro_noise_data gyro_noise];

        % y_noise=y_real;
        y_noise(4:6) = y_real(4:6) + gyro_noise;
        sign=randi([0 1]); 
        if sign==0 
            sign=-1; 
        end 
        y_noise(7:9) = css_noise(sun_pos_eci(:,(l-1)*10+1),x_real(1:4,(l-1)*10+1),xsat_eci(:,(l-1)*10+1),albedo(:,(l-1)*10+1),lambda);
        %y_noise(7:9) = css_noise(sun_pos_eci(:,(l-1)*10+1),x_real(1:4,(l-1)*10+1),xsat_eci(:,(l-1)*10+1),0,lambda);
        if eclipse((l-1)*10+1)~=0
            y_noise(7:9)=zeros(3,1);
        else
            y_noise(7:9)=y_noise(7:9)/norm(y_noise(7:9)); 
        end
        y_noise(1:3)=y_noise(1:3)/norm(y_noise(1:3)); 

        % y_noise=y_real;

        %% Wahba
        [q_wahba,~]=wahba(y_noise(7:9),y_noise(1:3),sun_pos_eci(:,(l-1)*10+1),mag_field_eci(:,(l-1)*10+1));
        ekf.theta(1:4) = q_wahba';


        %% Bias timer
        bias_init_counter = bias_init_counter + 1;
        % We add zeros for every timestep before the initialization is finished
        if (bias_init_counter < bias_wahba_loops + 1)        
            quat_pos(:,bias_init_counter) = q_wahba;    
            x_hat_data = [x_hat_data zeros(7,bias_wahba_loops)];  %Set measurements to 0 until bias has initialized
            
            % y_hat_data = [y_hat_data zeros(n_msr,10)]; 
            % P_data = [P_data zeros(10*10,10)];

            for i=1:10
                %% Propagate the system
                torq(1,1) = -(-(Const.I(3, 3) - Const.I(2, 2)) * x(6) * x(7))*(1 + 1e-2*randn())-kd*x(5);
                torq(2,1) = -(-(Const.I(1, 1) - Const.I(3, 3)) * x(5) * x(7))*(1 + 1e-2*randn())-kd*x(6);
                torq(3,1) = -(-(Const.I(2, 2) - Const.I(1, 1)) * x(5) * x(6))*(1 + 1e-2*randn())-kd*x(7);

                x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,rw_ang_momentum,[0;0;0])); 
                x_real(:,(l-1)*10+i+1)=x;
                t = t + dt;

            end
            continue
        end
        %% Bias calculation
        dqdt = zeros(4,bias_wahba_loops-1);  
        for i=2:bias_wahba_loops
            dqdt(:,i-1) = quat_pos(:,i)-quat_pos(:,i-1);
        end

        estimated_rate = zeros(3,bias_wahba_loops-1);

        % The angular rate is calculated using: angular_rate = 2 * dq/dt * q^-1
        for i=1:bias_wahba_loops-1
           temp = 2*quatProd(quatconj(quat_pos(:,i)'),dqdt(:,i));
           estimated_rate(:,i) = temp(2:4);
        end

        real_rate = y_noise(4:6);

        mean_omega = [0;0;0];        % Time averaging for noise reduction

        for i=1:3
           mean_omega(i) = mean(estimated_rate(i,:));  
        end

        initial_bias_estimate = real_rate-mean_omega;
        ekf.theta(5:7) = initial_bias_estimate; %Initialize angular velocity equal to gyroscope measurement 

        x_hat = ekf.theta;

        % Main continuous loop
        for k=l:n_steps
                % Artificial eclipse
    %     if k>2000
    %         for m=1:10
    %             eclipse(10*(k-1)+m)=2;
    %         end
    %     end

     q_ob_data(:,(k-1)*10+1) = quat_EB2OB(x(1:4), nodem(1,(k-1)*10+1),inclm(1,(k-1)*10+1),argpm(1,(k-1)*10+1),mm(1,(k-1)*10+1) );   
     for c=1:3
         
%% ================ Adaptive ======================================
%
%           if(entered_eclipse == false) && (eclipse((k-1)*10+c) ~= 0)
%               Q_eclipse =  eye(n_dim_error,n_dim_error);
%               Q_eclipse(4:6,4:6) = 0.5e-06*eye(3,3);
%               %Q_eclipse = Q_eclipse_load;
%               R_hat_coeff=10000*[.5e-3;.5e-3;.5e-3;4e-3;4e-3;4e-3;1e-3;1e-3;1e-3];
%               R_hat = R_hat_coeff.*eye(n_msr,n_msr);
%               ekf.setProcessNoiseCov(Q_eclipse); %Q variance matrix
%               ekf.setMeasureNoiseCov(R_hat); %R variance matrix
%               entered_eclipse = true;
%           end
%           
%           if(exited_eclipse == false) && (entered_eclipse == true) && (eclipse((k-1)*10+c) == 0)
%               Q_outside =  0.5e-05*eye(n_dim_error,n_dim_error);
%               Q_outside(4:6,4:6) = 0.5e-07*eye(3,3);
%               %Q_outside = Q_eclipse_load;
%               %Q_outside = Q;
%               ekf.setProcessNoiseCov(Q_outside);
%               ekf.setMeasureNoiseCov(Param.R_hat); %R variance matrix
%               exited_eclipse = true;
%               entered_eclipse = false;
%           end
%% ============ IdealQ ============================================        
        if (eclipse((k-1)*10+c))
            % Variances
            Q = Q_eclipse_load; % Variance of the process noise w[k]

            % R Variances used in EKF
            % R_hat_coeff=[1e-3;1e-3;1e-3;8e-3;8e-3;8e-3;5e-3;5e-3;5e-3];
            R_hat_coeff=10000*[.5e-3;.5e-3;.5e-3;4e-3;4e-3;4e-3;1e-3;1e-3;1e-3];
            R_hat = R_hat_coeff.*eye(n_msr,n_msr);
            ekf.setProcessNoiseCov(Q); %Q variance matrix
             ekf.setMeasureNoiseCov(R_hat); %R variance matrix
        else
            ekf.setProcessNoiseCov(Param.Q); %Q variance matrix
            ekf.setMeasureNoiseCov(Param.R_hat); %R variance matrix
        end
       
        
        y_real = real_model.msrFun(x_real(:,(k-1)*10+c),msrCookieFinal(mag_field_eci(:,(k-1)*10+c),...
            sun_pos_eci(:,(k-1)*10+c),eclipse((k-1)*10+c),[0;0;0]));
        Bbody_data(:,(k-1)*10+c) = y_real(1:3)*norm(mag_field_orbit(:,(k-1)*10+c)*10^(-9));
        y_noise = y_real + sqrt(R)*randn(size(y_real));
        [gyro_noise,real_bias] = gyro_noise_func(real_bias,dt,sigma_u,sigma_v);
        bias_data = [bias_data real_bias];
        gyro_noise_data = [gyro_noise_data gyro_noise];
        y_noise(4:6) = y_real(4:6) + gyro_noise;
             
        y_noise(7:9) = css_noise(sun_pos_eci(:,(k-1)*10+c),x_real(1:4,(k-1)*10+c),xsat_eci(:,(k-1)*10+c),albedo(:,(k-1)*10+c),lambda);
        %y_noise(7:9) = css_noise(sun_pos_eci(:,(k-1)*10+c),x_real(1:4,(k-1)*10+c),xsat_eci(:,(k-1)*10+c),0,lambda);      
        
        if eclipse((k-1)*10+c)~=0
            y_noise(7:9)=zeros(3,1);
        else
            y_noise(7:9)=y_noise(7:9)/norm(y_noise(7:9)); 
        end
        y_noise(1:3)=y_noise(1:3)/norm(y_noise(1:3)); 
        
        gyro = y_noise(4:6);
        ekf.correct(y_noise, msrCookieFinalExtended(mag_field_eci(:,(k-1)*10+c),...
            sun_pos_eci(:,(k-1)*10+c),eclipse((k-1)*10+c),gyro,xsat_eci(:,(k-1)*10+c),albedo_inaccurate(:,(k-1)*10+c),lambda));
            %sun_pos_eci(:,(k-1)*10+c),eclipse((k-1)*10+c),gyro,xsat_eci(:,(k-1)*10+c),0,lambda));
        
        x_hat = ekf.theta;
        x_hat(1:4) = x_hat(1:4) / norm(x_hat(1:4));
        x_hat_data = [x_hat_data x_hat];

        %% Propagate the system
        torq=zeros(3,1);
        
        [~, ~, ~, AngVel_rw_radps(3,1), acceleration_rw(2,1), T_rw_total] = rw_model(0, AngVel_rw_radps(3,1)); 
        torq(3,1) = T_rw_total;
        
        [T_dist, residual,~,ad,r,sp,g] = disturbances_pd(q_ob_data(:,(k-1)*10+c), sun_pos_orbit(:,(k-1)*10+c), mag_field_orbit(:,(k-1)*10+c)*10^(-9), disturbancesEnabled);
        rm = [rm residual];
        tau_dist(:,(k-1)*10+c+1) = T_dist;
        torq = torq + T_dist;
        tau_ad(:,(k-1)*10+c+1) = ad;
        tau_rm(:,(k-1)*10+c+1) = r;
        tau_sp(:,(k-1)*10+c+1) = sp;
        tau_g(:,(k-1)*10+c+1) = g;
        tau_dist(:,(k-1)*10+c+1) = T_dist;
        torq = torq + T_dist;
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,rw_ang_momentum,[0;0;0])); 
        x_real(:,(k-1)*10+c+1)=x;
        q_ob_data(:,(k-1)*10+c+1) = quat_EB2OB(x(1:4), nodem(1,(k-1)*10+c),inclm(1,(k-1)*10+c),argpm(1,(k-1)*10+c),mm(1,(k-1)*10+c) );
        
        AngVel_rw_rpm(3,1) = AngVel_rw_radps(3,1)*30/pi; 
        rw_ang_vel_rpm(1,(k-1)*10+c+1) = AngVel_rw_rpm(3,1);  
        %rw_accel(1,(k-1)*10+c+1) = acceleration_rw(2,1);
        
    %     if norm(q_prev - q_ob_data(:,(k-1)*10+c+1))^2 > norm(q_prev + q_ob_data(:,(k-1)*10+c+1))^2
    %         q_ob_data(:,(k-1)*10+c+1) = -q_ob_data(:,(k-1)*10+c+1);
    %     end
%          if q_ob_data(1,(k-1)*10+c+1)<0
%             q_ob_data(:,(k-1)*10+c+1) = -q_ob_data(:,(k-1)*10+c+1);
%          end  
    %    q_prev = q_ob_data(:,(k-1)*10+c+1);
        % Predict the states at next time step, k+1. This updates the State and
        % StateCovariance properties of the filter to contain x[k+1|k] and
        % P[k+1|k]. These will be utilized by the filter at the next time step.
        
        gyro = y_noise(4:6);
        ekf.predict(stateTransCookieFinalNominal(torq,rw_ang_momentum,gyro),dt);
     end 
    %  
    %  q_ob_hat_prev =quat_EB2OB(x_hat(1:4), nodem(1,(k-1)*10+c),inclm(1,(k-1)*10+c),argpm(1,(k-1)*10+c),mm(1,(k-1)*10+c) );
    % 
    % if q_ob_hat_prev(1) < 0
    %    q_ob_hat_prev = -q_ob_hat_prev;
    % end
     for c=4:10   
        
%% ================ Adaptive ======================================
%
%           if(entered_eclipse == false) && (eclipse((k-1)*10+c) ~= 0)
%               Q_eclipse =  eye(n_dim_error,n_dim_error);
%               Q_eclipse(4:6,4:6) = 0.5e-06*eye(3,3);
%               %Q_eclipse = Q_eclipse_load;
%               R_hat_coeff=10000*[.5e-3;.5e-3;.5e-3;4e-3;4e-3;4e-3;1e-3;1e-3;1e-3];
%               R_hat = R_hat_coeff.*eye(n_msr,n_msr);
%               ekf.setProcessNoiseCov(Q_eclipse); %Q variance matrix
%               ekf.setMeasureNoiseCov(R_hat); %R variance matrix
%               entered_eclipse = true;
%           end
%           
%           if(exited_eclipse == false) && (entered_eclipse == true) && (eclipse((k-1)*10+c) == 0)
%               Q_outside =  0.5e-05*eye(n_dim_error,n_dim_error);
%               Q_outside(4:6,4:6) = 0.5e-07*eye(3,3);
%               %Q_outside = Q_eclipse_load;
%               %Q_outside = Q;
%               ekf.setProcessNoiseCov(Q_outside);
%               ekf.setMeasureNoiseCov(Param.R_hat); %R variance matrix
%               exited_eclipse = true;
%               entered_eclipse = false;
%           end
%% ============ IdealQ ============================================        
        if (eclipse((k-1)*10+c))
            % Variances
            Q = Q_eclipse_load; % Variance of the process noise w[k]

            % R Variances used in EKF
            % R_hat_coeff=[1e-3;1e-3;1e-3;8e-3;8e-3;8e-3;5e-3;5e-3;5e-3];
            R_hat_coeff=10000*[.5e-3;.5e-3;.5e-3;4e-3;4e-3;4e-3;1e-3;1e-3;1e-3];
            R_hat = R_hat_coeff.*eye(n_msr,n_msr);
            ekf.setProcessNoiseCov(Q); %Q variance matrix
             ekf.setMeasureNoiseCov(R_hat); %R variance matrix
        else
            ekf.setProcessNoiseCov(Param.Q); %Q variance matrix
            ekf.setMeasureNoiseCov(Param.R_hat); %R variance matrix
        end

        y_real = real_model.msrFun(x_real(:,(k-1)*10+c),msrCookieFinal(mag_field_eci(:,(k-1)*10+c),...
            sun_pos_eci(:,(k-1)*10+c),eclipse((k-1)*10+c),[0;0;0]));
        Bbody_data(:,(k-1)*10+c) = y_real(1:3)*norm(mag_field_orbit(:,(k-1)*10+c)*10^(-9));
        [gyro_noise,real_bias] = gyro_noise_func(real_bias,dt,sigma_u,sigma_v);
        bias_data = [bias_data real_bias];
        gyro_noise_data = [gyro_noise_data gyro_noise];
        y_noise(4:6) = y_real(4:6) + gyro_noise; 
        
        x_hat = ekf.theta;
        x_hat(1:4) = x_hat(1:4) / norm(x_hat(1:4));
        x_hat_data = [x_hat_data x_hat];
        q_ob_hat = quat_EB2OB(x_hat(1:4), nodem(1,(k-1)*10+c-1),inclm(1,(k-1)*10+c-1),argpm(1,(k-1)*10+c-1),mm(1,(k-1)*10+c-1) );

    %     if norm(q_ob_hat_prev - q_ob_hat)^2 > norm(q_ob_hat_prev + q_ob_hat)^2
    %        q_ob_hat = -q_ob_hat;
    %     end
%         if q_ob_hat(1) < 0
%             q_ob_hat = -q_ob_hat;
%         end
    %     
        q_ob_hat_prev = q_ob_hat;

        
        %% PD function
        acceleration_rw(1,1) = acceleration_rw(2,1);
        acceleration_rw(2,1) = acceleration_rw(3,1);
        AngVel_rw_radps(1,1) = AngVel_rw_radps(2,1);
        AngVel_rw_radps(2,1) = AngVel_rw_radps(3,1);
        AngVel_rw_rpm(1,1) = AngVel_rw_rpm(2,1);
        AngVel_rw_rpm(2,1) = AngVel_rw_rpm(3,1);
        
        [torq, T_rw, T_magnetic_effective, V_rw, I_rw, P_thermal_rw, AngVel_rw_rpm_next, AngVel_rw_radps_next,...
                acceleration_rw_cur, rw_ang_momentum, init_AngVel_dz, init_accel_dz, V_mtq, I_mtq, P_thermal_mtq, ...
                    timeflag_dz,M] = ...
                        PD(Kp_gain, Kd_gain, q_desired ,q_ob_hat, y_noise(4:6)-ekf.theta(5:7) , y_noise(1:3)*norm(mag_field_orbit(:,(k-1)*10+c)*10^(-9)) , eclipse((k-1)*10+c), ...
                            Const.mtq_max1, Const.mtq_max2, Const.mtq_max3, Const.lim_dz, AngVel_rw_radps(2,1), AngVel_rw_rpm(2,1), ...
                                acceleration_rw(1,1), init_AngVel_dz, init_accel_dz, timeflag_dz,Const.rw_max_torque,y_real(1:3)*norm(mag_field_orbit(:,(k-1)*10+c)*10^(-9)), l, sat_llh(3,(k-1)*10+c));
        tau_rw(1, (k-1)*10+c+1) = T_rw(3);
        tau_mtq(:, (k-1)*10+c+1) = T_magnetic_effective;
        rw_ang_vel_rpm(1,(k-1)*10+c+1) = AngVel_rw_rpm(3,1); 
        acceleration_rw(2,1) = acceleration_rw_cur;
        AngVel_rw_rpm(3,1) = AngVel_rw_rpm_next;
        AngVel_rw_radps(3,1) = AngVel_rw_radps_next;
        M_data(:,(k-1)*10+c+1) = M;
        %%                    
        [T_dist, residual,~,ad,r,sp,g] = disturbances_pd(q_ob_data(:,(k-1)*10+c), sun_pos_orbit(:,(k-1)*10+c), mag_field_orbit(:,(k-1)*10+c)*10^(-9), disturbancesEnabled);
        rm = [rm residual];
        tau_dist(:,(k-1)*10+c+1) = T_dist;
        torq = torq + T_dist;
        tau_ad(:,(k-1)*10+c+1) = ad;
        tau_rm(:,(k-1)*10+c+1) = r;
        tau_sp(:,(k-1)*10+c+1) = sp;
        tau_g(:,(k-1)*10+c+1) = g;
        tau_dist(:,(k-1)*10+c+1) = T_dist;
        torq = torq + T_dist;
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,rw_ang_momentum,[0;0;0])); 
        x_real(:,(k-1)*10+c+1)=x;
        q_ob_data(:,(k-1)*10+c+1) = quat_EB2OB(x(1:4), nodem(1,(k-1)*10+c),inclm(1,(k-1)*10+c),argpm(1,(k-1)*10+c),mm(1,(k-1)*10+c) );



    %     if norm(q_prev - q_ob_data(:,(k-1)*10+c+1))^2 > norm(q_prev + q_ob_data(:,(k-1)*10+c+1))^2
    %         q_ob_data(:,(k-1)*10+c+1) = -q_ob_data(:,(k-1)*10+c+1);
    %     end
%         if q_ob_data(1,(k-1)*10+c+1)<0
%             q_ob_data(:,(k-1)*10+c+1) = -q_ob_data(:,(k-1)*10+c+1);
%         end  
    %     q_prev = q_ob_data(:,(k-1)*10+c+1);
        gyro = y_noise(4:6);
        ekf.predict(stateTransCookieFinalNominal(torq,rw_ang_momentum,gyro),dt);
     end
        end
        break;
    else
        for i=1:10
            %% Propagate the system
            torq=zeros(3,1);

            x = real_model.stateTransFun(x, stateTransCookieFinal(torq)); 
            x_real(:,(l-1)*10+i+1)=x;
            t = t + dt;

        end
        x_hat_data = [x_hat_data zeros(10,10)];
        bias_data = [bias_data zeros(3,1)];
        gyro_noise_data = [gyro_noise_data zeros(3,1)];
    end

end

%% =============================== Errors and plots =============================================== %%

%% Calculation and plotting of performance error
x_real_euler_perf = quat2eul(q_ob_data(1:4,1:length(q_ob_data))'); 
x_real_euler_perf = rad2deg(x_real_euler_perf');

instant_error_perform = x_real_euler_perf';


figure();
for i=1:3
    subplot(3,1,i);
    hold on;
    plot(Time(21:length(instant_error_perform)), instant_error_perform(21:length(instant_error_perform), i), 'LineWidth',1.5, 'Color','blue');
    if (i==1), title('Absolute Performance Errors', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    if (i==2), ylabel('Y-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    if (i==3), ylabel('Z-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

%% Calulation and plotting of knowledge error


     
x_hat_euler_know = zeros(length(x_hat_data), 6);
instant_error_know = zeros(length(x_hat_data), 6);

x_real_euler_know = quat2eul(x_real(1:4,1:length(x_hat_data))');
x_real_euler_know = rad2deg(x_real_euler_know');
x_hat_euler_know(:, 1:3) = quat2eul(x_hat_data(1:4,:)');
x_hat_euler_know(:, 1:3) = (rad2deg(x_hat_euler_know(:, 1:3)'))';

instant_error_know(:, 1:3) = x_hat_euler_know(:, 1:3) - x_real_euler_know';
instant_error_know(:, 4:6) = x_hat_data(5:7, 1:length(x_hat_data))' - bias_data';

for i=1:3
    for k=1:length(instant_error_know)
        if instant_error_know(k, i) > 180
            instant_error_know(k, i) = instant_error_know(k, i) - 360;
        elseif instant_error_know(k, i) < -180
            instant_error_know(k, i) = instant_error_know(k, i) + 360;
        end
    end
end

figure();
for i=1:6
    subplot(6,1,i);
    hold on;
    plot(Time(21:length(instant_error_know)), instant_error_know(21:length(instant_error_know), i), 'LineWidth',1.5, 'Color','blue');
    if (i==1), title('Absolute Knowledge Errors', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    if (i==2), ylabel('Y-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    if (i==3), ylabel('Z-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    if (i==4) ylabel(['$\omega_1 [rad/sec]$'], 'interpreter','latex', 'fontsize',14); end
    if (i==5) ylabel(['$\omega_2 [rad/sec]$'], 'interpreter','latex', 'fontsize',14); end
    if (i==6) ylabel(['$\omega_3 [rad/sec]$'], 'interpreter','latex', 'fontsize',14); end
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

n_dim = size(x_real,1);
 %figure('Position',[500 0 1420 1080]);
figure();
for i=1:n_dim
    subplot(n_dim,1,i);
    hold on;
    if i<5
        plot(Time,x_real(i,1:length(Time)), 'LineWidth',2.0, 'Color','blue');       
    else
        plot(Time(1:length(bias_data(1,:))),bias_data(i-4,1:length(bias_data(1,:))))
    end
    
    plot(Time(1:length(x_hat_data(i,:))),x_hat_data(i,:), 'LineWidth',2.0, 'Color','magenta');               
    if (i==1),legend({['$x_' num2str(i) '$'],['$\hat{x}_' num2str(i) '$']}, 'interpreter','latex', 'fontsize',15);end
    ylabel(['$x_' num2str(i) '$'], 'interpreter','latex', 'fontsize',14);
    if (i==1), title('MEKF estimation results', 'interpreter','latex', 'fontsize',17);end
%     xlim([3 n_steps]);
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

% % figure('Position',[1000 0 1420 1080]);
for i=1:length(Time)
    if(q_ob_data(1,i)<0)
    q_ob_data(:,i) = -q_ob_data(:,i);
    end
end

figure();
for i=1:4
    subplot(4,1,i);
    if (i==1), title('Quaternion', 'interpreter','latex', 'fontsize',17); end
    hold on;
    plot(Time,q_ob_data(i,1:length(Time)), 'LineWidth',2.0, 'Color','blue');
    ylabel(['$q_{ob' num2str(i) '}$'], 'interpreter','latex', 'fontsize',14);
    xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
%     xlim([3 n_steps]);
    hold off;
    grid on;
end

%% Eclipse plot

    figure()
    plot(1:length(eclipse),eclipse, 'LineWidth',2.0, 'Color','blue');
    title('Eclipse', 'interpreter','latex', 'fontsize',17)
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    ylabel(['Eclipse'], 'interpreter','latex', 'fontsize',14);
    grid on;
    if (i==1), title('Umbral, Penumbral or no Eclipse', 'interpreter','latex', 'fontsize',17);end

%%
% x_err_data(1:4,:) = x_real(1:4,1:length(x_hat_data))-x_hat_data(1:4,:);
% x_err_data(5:7,:) = bias_data - x_hat_data(5:7,:);
% 
% figure();
% for i=1:n_dim
%     subplot(n_dim,1,i);
%     plot(Time(1:length(x_err_data(i,:))),x_err_data(i,:), 'LineWidth',2.0, 'Color','blue');  % Error for the first state
%     xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
%     ylabel(['$\tilde{x}_' num2str(i) '$'], 'interpreter','latex', 'fontsize',14);
%     if (i==1), legend({'State estimate','$\pm \sigma$'}, 'interpreter','latex', 'fontsize',11);end
%     if (i==1), title('State estimation errors', 'interpreter','latex', 'fontsize',11); end
%     xlim([3 n_steps]);
%     hold off;
% end
%

%%
    figure();
    for i=1:3
        subplot(3,1,i);
        plot(Time,x_real(4+i,1:length(Time)),'LineWidth',2.0, 'Color','blue');
        xlabel('Time [s]', 'interpreter','latex', 'fontsize',12);
        if (i==1) ylabel(['\omega1 [rad/sec]'], 'interpreter','latex', 'fontsize',14); end
        if (i==2) ylabel(['\omega2 [rad/sec]'], 'interpreter','latex', 'fontsize',14); end
        if (i==3) ylabel(['\omega3 [rad/sec]'], 'interpreter','latex', 'fontsize',14); end
        ylabel(['$\omega_' num2str(i) '$' '[rad/sec]'], 'interpreter','latex', 'fontsize',14);
        if (i==1), legend('Angular Velocity');end
        if (i==1), title('Angular Velocities', 'interpreter','latex', 'fontsize',17); end
        grid on;
    end
%%
% 
% figure();
% for i=1:3
%     subplot(3,1,i);
%     plot(Time(1:length(x_hat_data(1,:))), x_hat_data(4+i,:) - gyro_noise_data(i,:),'LineWidth',2.0, 'Color','blue');
%     xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
%     ylabel(['$\omega_' num2str(i) '$'], 'interpreter','latex', 'fontsize',14);
%     if (i==1), legend('Angular Velocity estimation error');end
%     if (i==1), title('Angular Velocity estimation error'); end
% end
% 
% eul_diff =zeros(length(Time),3);
% euler_hat=quat2eul(q_ob_data(1:4,:)')';
% eul_diff=(euler_hat)*180/pi;
% figure();
% for i=1:3
%     subplot(3,1,i);
%     hold on;
%     plot(Time,eul_diff(i,1:length(Time)), 'LineWidth',2.0, 'Color','blue');
%     ylabel(['$euler_{' num2str(i) '}$'], 'interpreter','latex', 'fontsize',17);
%     xlim([3 n_steps]);
%     hold off;
% end
%%  Plotting the produced Torques

      total_torques = tau_mtq;
      total_torques(3,:) = total_torques(3,:) + tau_rw;
 
      figure()
      subplot(3,3,1)
      plot(1:length(Time),tau_mtq(1,1:length(Time)))
      title('Magnetic Torques', 'interpreter','latex', 'fontsize',17)
      ylabel('Torque-x [Nm]', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;
      subplot(3,3,4)
      plot(1:length(Time),tau_mtq(2,1:length(Time)))
      ylabel('Torque-y [Nm]', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;
      subplot(3,3,7)
      plot(1:length(Time),tau_mtq(3,1:length(Time)))
      ylabel('Torque-z [Nm]', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;
      
      subplot(3,3,3)
      plot(1:length(Time),total_torques(1,1:length(Time)))
      ylabel('Torques - x', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      title('Total Torques', 'interpreter','latex', 'fontsize',17)
      grid on;
      
      subplot(3,3,6)
      plot(1:length(Time),total_torques(2,1:length(Time)))
      ylabel('Torques - y', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;      
      
      subplot(3,3,9)
      plot(1:length(Time),total_torques(3,1:length(Time)))
      ylabel('Torques - z', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;

      subplot(3,3,8)
      plot(1:length(Time),tau_rw(1, 1:length(Time)))
      title('Reaction Wheel Torque-z', 'interpreter','latex', 'fontsize',17)
      ylabel('Torque [Nm]', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;

      subplot(3,3,2)
      plot(Time, rw_ang_vel_rpm(1,1:length(Time)),'LineWidth',1.5, 'Color','blue'); 
      title('Angular velocity of RW', 'interpreter','latex', 'fontsize',17); 
      ylabel('Angular Velocity [rpm]', 'interpreter','latex', 'fontsize',14); 
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12); 
      grid on;
      
%% 
mtq_max_vector = [Const.mtq_max1 Const.mtq_max1 Const.mtq_max2]';
mtq_max_vector = mtq_max_vector * ones(1,length(Bbody_data));
max_mtq_torque = abs(cross(mtq_max_vector,Bbody_data));
figure()
for i=1:3
    subplot(3,3,i);
    hold on;
    plot(Time(21:length(rm)), max_mtq_torque(i, 21:length(rm)), 'LineWidth',1.5, 'Color','blue');
    if (i==2), title('Maximum MTQ torque', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis [deg]'); end
    if (i==2), ylabel('Y-axis [deg]'); end
    if (i==3), ylabel('Z-axis [deg]'); end
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

for i=1:3
    subplot(3,3,i+3);
    hold on;
    plot(Time(21:length(rm)), total_torques(i,21:length(rm)), 'LineWidth',1.5, 'Color','blue');
    if (i==2), title('Total actuator torque', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis [deg]'); end
    if (i==2), ylabel('Y-axis [deg]'); end
    if (i==3), ylabel('Z-axis [deg]'); end
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end
for i=1:3
    subplot(3,3,i+6);
    hold on;
    plot(Time(22:length(rm)),tau_dist(i,22:length(rm)), 'LineWidth',1.5, 'Color','blue')
    if (i==2), title('Total disturbance torque', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis [deg]'); end
    if (i==2), ylabel('Y-axis [deg]'); end
    if (i==3), ylabel('Z-axis [deg]'); end
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

%Max actuator Torque - Actuator Torque - Total Disturbances
      
      
      figure()
      subplot(3,1,1)
      plot(1:length(Time),M_data(1, 1:length(Time)))
      title('Induced dipole - x', 'interpreter','latex', 'fontsize',17)
      ylabel('Torque [Nm]', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;
      subplot(3,1,2)
      plot(1:length(Time),M_data(2, 1:length(Time)))
      title('Induced dipole - y', 'interpreter','latex', 'fontsize',17)
      ylabel('Torque [Nm]', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;
      subplot(3,1,3)
      plot(1:length(Time),M_data(3, 1:length(Time)))
      title('Induced dipole - z', 'interpreter','latex', 'fontsize',17)
      ylabel('Torque [Nm]', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;
%% 
%      figure() 
%      plot(Time, rw_ang_vel_rpm(1,1:length(Time)),'LineWidth',1.5, 'Color','blue'); 
%      title('Angular velocity of RW', 'interpreter','latex', 'fontsize',17); 
%      ylabel('Angular Velocity [rpm]', 'interpreter','latex', 'fontsize',14); 
%      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12); 
%      grid on;
 
 %%  Plotting the Disturbances
 
      figure()
      subplot(3,1,1)
      plot(1:length(Time),tau_dist(1,1:length(Time)), 'LineWidth',1.5, 'Color','blue')
      title('Disturbances', 'interpreter','latex', 'fontsize',17)
      ylabel('Disturbances-x [Nm]', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;
      subplot(3,1,2)
      plot(1:length(Time),tau_dist(2,1:length(Time)), 'LineWidth',1.5, 'Color','blue')
      ylabel('Disturbances-y [Nm]', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;
      subplot(3,1,3)
      plot(1:length(Time),tau_dist(3,1:length(Time)), 'LineWidth',1.5, 'Color','blue')
      ylabel('Disturbances-z [Nm]', 'interpreter','latex', 'fontsize',14)
      xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
      grid on;
      
%% Calculation and plotting of Mean Performance Error

mean_error_perf = zeros(3, 11);
for i = 1:3
    for j = 1:11*orbits
    mean_error_perf(i, j) = mean(instant_error_perform(21+((500/dt)*(j-1)):21+((500/dt)*j), i));
    end
end

mean_error_perf_matrix = zeros(length(x_hat_data),3);
for j = 1:11*orbits
    for i = 21+((500/dt)*(j-1)):21+((500/dt)*j)
        mean_error_perf_matrix(i, :) = mean_error_perf(:,j);
    end
end

figure();
for i=1:3
    subplot(3,1,i);
    hold on;
    plot(Time(1:length(mean_error_perf_matrix)), mean_error_perf_matrix(1:length(mean_error_perf_matrix), i), 'LineWidth',1.5, 'Color','blue');
    if (i==1), title('Mean Performance Errors', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    if (i==2), ylabel('Y-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    if (i==3), ylabel('Z-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    xlabel('Time [s]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

 %% Calculation and plotting of Mean Knowledge Error

mean_error_know = zeros(6, 11);
for i = 1:6
    for j = 1:11*orbits
    mean_error_know(i, j) = mean(instant_error_know(21+((500/dt)*(j-1)):21+((500/dt)*j), i));
    end
end

mean_error_know_matrix = zeros(length(x_hat_data),6);
for j = 1:11*orbits
    for i = 21+((500/dt)*(j-1)):21+((500/dt)*j)
        mean_error_know_matrix(i, :) = mean_error_know(:,j);
    end
end

figure();
for i=1:3
    subplot(3,1,i);
    hold on;
    plot(Time(1:length(mean_error_know_matrix)), mean_error_know_matrix(1:length(mean_error_know_matrix), i), 'LineWidth',1.5, 'Color','blue');
    if (i==1), title('Mean Knowledge Errors', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    if (i==2), ylabel('Y-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    if (i==3), ylabel('Z-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    xlabel('Time [s]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end 


%% Calculation and plotting of Relative Performance Error
relative_error_perf = zeros(length(mean_error_perf_matrix(:,i)),3);
for i=1:3
relative_error_perf(:,i) = instant_error_perform(1:length(mean_error_perf_matrix(:,i)),i) - mean_error_perf_matrix(:,i);
end

figure();
for i=1:3
    subplot(3,1,i);
    hold on;
    plot(Time(21:length(relative_error_perf)), relative_error_perf(21:length(relative_error_perf), i), 'LineWidth',1.5, 'Color','blue');
    if (i==1), title('Relative Performance Errors', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis [deg]'); end
    if (i==2), ylabel('Y-axis [deg]'); end
    if (i==3), ylabel('Z-axis [deg]'); end
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

%% Calculation and plotting of Relative Knowledge Error
relative_error_know = zeros(length(instant_error_know),6);
for i=1:6
relative_error_know(:,i) = instant_error_know(:,i) - mean_error_know_matrix(1:length(instant_error_know(:,i)),i);
end

figure();
for i=1:6
    subplot(6,1,i);
    hold on;
    plot(Time(21:length(instant_error_know)), relative_error_know(21:length(relative_error_know), i), 'LineWidth',1.5, 'Color','blue');
    ylabel(['$\tilde{x}_' num2str(i) ' [rad/sec]$'], 'interpreter','latex', 'fontsize',14);
    if (i==1), title('Relative Knowledge Errors', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    if (i==2), ylabel('Y-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    if (i==3), ylabel('Z-axis [deg]', 'interpreter','latex', 'fontsize',14); end
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

figure();
for i=1:3
    subplot(3,1,i);
    hold on;
    plot(Time(21:length(rm)), rm(i, 21:length(rm)), 'LineWidth',1.5, 'Color','blue');
    if (i==1), title('Residual Magnetic Dipole', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis [deg]'); end
    if (i==2), ylabel('Y-axis [deg]'); end
    if (i==3), ylabel('Z-axis [deg]'); end
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

%% Plotting disturbances
 
    
     T_disturbances = tau_rm + tau_ad + tau_sp + tau_g;
     
     subplot(3,5,1)
     plot(1:length(Time),T_disturbances(1,1:length(Time)))
     ylabel('Torque [Nm]')
     xlabel('Time [s]');
     title('Disturbances-x')
     grid on;
     subplot(3,5,6)
     plot(1:1:length(Time),T_disturbances(2,1:length(Time)))
     ylabel('Torque [Nm]')
     title('Disturbances-y')
     xlabel('Time [s]');
     grid on;
     subplot(3,5,11)
     plot(1:1:length(Time),T_disturbances(3,1:length(Time)))
     title('Disturbances-z')
     ylabel('Torque [Nm]')
     xlabel('Time [s]');
     grid on;
     
     %% Plotting rm
 

     subplot(3,5,2)
     plot(1:length(Time),tau_rm(1,1:length(Time)))
     ylabel('Torque [Am^2]')
     xlabel('Time [s]');
     title('MRD-x')
     grid on;
     subplot(3,5,7)
     plot(1:length(Time),tau_rm(2,1:plotter_step:end))
     ylabel('Torque [Am^2]')
     title('MRD-y')
     xlabel('Time [s]');
     grid on;
     subplot(3,5,12)
     plot(1:plotter_step:reps,tau_rm(3,1:plotter_step:end))
     title('MRD-z')
     ylabel('Torque [Am^2]')
     xlabel('Time [s]');
     grid on;
     
     %% Plotting ad
 
    
     subplot(3,5,3)
     plot(1:plotter_step:reps,tau_ad(1,1:plotter_step:end))
     ylabel('Torque [Nm]')
     xlabel('Time [s]');
     title('Aerodynamic-x')
     grid on;
     subplot(3,5,8)
     plot(1:plotter_step:reps,tau_ad(2,1:plotter_step:end))
     ylabel('Torque [Nm]')
     title('Aerodynamic-y')
     xlabel('Time [s]');
     grid on;
     subplot(3,5,13)
     plot(1:plotter_step:reps,tau_ad(3,1:plotter_step:end))
     title('Aerodynamic-z')
     ylabel('Torque [Nm]')
     xlabel('Time [s]');
     grid on;
     
     %% Plotting g
 
 
     subplot(3,5,4)
     plot(1:plotter_step:reps,tau_g(1,1:plotter_step:end))
     ylabel('Torque [Nm]')
     xlabel('Time [s]');
     title('Gravitational-x')
     grid on;
     subplot(3,5,9)
     plot(1:plotter_step:reps,tau_g(2,1:plotter_step:end))
     ylabel('Torque [Nm]')
     title('Gravitational-y')
     xlabel('Time [s]');
     grid on;
     subplot(3,5,14)
     plot(1:plotter_step:reps,tau_g(3,1:plotter_step:end))
     title('Gravitational-z')
     ylabel('Torque [Nm]')
     xlabel('Time [s]');
     grid on;
     
     %% Plotting sp
 
 
     subplot(3,5,5)
     plot(1:plotter_step:reps,tau_sp(1,1:plotter_step:end))
     ylabel('Torque [Nm]')
     xlabel('Time [s]');
     title('Solar Pressure-x')
     grid on;
     subplot(3,5,10)
     plot(1:plotter_step:reps,tau_sp(2,1:plotter_step:end))
     ylabel('Torque [Nm]')
     title('Solar Pressure-y')
     xlabel('Time [s]');
     grid on;
     subplot(3,5,15)
     plot(1:plotter_step:reps,tau_sp(3,1:plotter_step:end))
     title('Solar Pressure-z')
     ylabel('Torque [Nm]')
     xlabel('Time [s]');
     grid on;






end