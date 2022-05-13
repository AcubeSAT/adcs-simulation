function measurement_fault_detection
    close all;
    clc;
    
    Kp_gain= 1e-05*diag([20 150 120]); % 500km 
    Kd_gain= 1e-04*diag([75 100 75]);


%% Initialize Parameters Script

    Const=constants();
    Param = setParamsFinal_Nominal(Const.I);
    dt = Param.dt;
    orbits = Param.orbits;
    tf = Param.tf;
    q_desired = Param.q_desired;
    x0 = Param.x0;
    x0_hat = Param.x0_hat;
    real_model = Param.real_model;
    model = Param.model;

    albedo = Param.albedo;
    albedo_inaccurate = Param.albedo_inaccurate;
    disturbancesEnabled = Param.disturbancesEnabled;

    eclipse = Param.eclipse;
    mag_field_eci = Param.mag_field_eci;
    mag_field_orbit = Param.mag_field_orbit;
    sun_pos_eci = Param.sun_pos_eci;
    sun_pos_orbit = Param.sun_pos_orbit;
    xsat_eci = Param.xsat_eci;
    argpm = Param.argpm;
    nodem = Param.nodem;
    inclm = Param.inclm;
    mm = Param.mm;
    init_bias = Param.init_bias;
    Q = Param.Q;
    R = Param.R;
    R_hat = Param.R_hat;
    sigma_u = Param.sigma_u;
    sigma_v = Param.sigma_v;
    P0 = Param.P0;
    number_of_measurements = Param.number_of_measurements;
    use_analytic_jacob = Param.use_analytic_jacob;


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

%% Construct the MEKF

    Q_struct = load('idealQ.mat', 'IDEALQ');
    Q_eclipse_load = Q_struct.IDEALQ;
    n_params = length(x0_hat);                                                         % Init number of parameters
    mekf = MEKF(n_params, number_of_measurements, @model.stateTransFun, @model.msrFun); % Init EKF
    mekf.global_state = x0_hat;                                                        % Init state estimation
    mekf.P = P0;                                                                       % Init Covariance matrix
    mekf.setProcessNoiseCov(Q);                                                        % Q variance matrix
    mekf.setMeasureNoiseCov(R_hat);                                                    % R variance matrix
    mekf.setFadingMemoryCoeff(1.00);                                                   % This parameter defined the "memory" of the filter, 1 meaning regular
    mekf.setPartDerivStep(0.001);                                                      % Arithmetical jacobian step
    if (use_analytic_jacob)                                                            % If analytical jacobian is used, define the functions
        mekf.setStateTransFunJacob(@model.stateTransFunJacob);
        mekf.setMsrFunJacob(@model.msrFunJacob);
    end



%% Initialize simulation parameters

    Time = 0:dt:tf;
    x_real = zeros(7,length(Time));             % Real state
    q_ob_data = zeros(4,length(Time));
    x = x0(1:7);
    x_real(:,1) = x0(1:7);
    t = 0;
    N_Timesteps = 10;                           % Number of timesteps per cycle
    number_of_cycles = floor(length(Time)/N_Timesteps);  % Number of cycles
    timeflag_dz = 0;
    init_AngVel_dz = 0;
    init_accel_dz = 0;
    rw_ang_momentum=0;
    rw_ang_vel_rpm = zeros(1, length(Time));    % RW angular velocity in rpm
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
    AngVel_rw_radps = zeros(3, 1);              % 1 = old, 2 = cur, 3 = next
    AngVel_rw_rpm = zeros(3, 1);
    acceleration_rw = zeros(3, 1);
    x_hat_data = zeros(7,length(Time));
    bias_data = zeros(3,length(Time));
    gyro_noise_data = zeros(3,length(Time));
    Bbody_data = zeros(3,length(Time));
    bdot_activation_matrix = zeros(2, length(Time));
    threshold_times = 0; 
    threshold_exceptions = 0;
    
%% Next we initialize the bias estimation by solving Wahba's problem n times. 

    bias_init_counter = 0;                      % how many Wahba's have we solved for bias init
    bias_wahba_loops = 2;                       % Total times to be solved
    quat_pos = zeros(4,bias_wahba_loops);       % Wahba results are stored here
    real_bias=init_bias;


    for cycle_index = 1:bias_wahba_loops
        current_timestep = (cycle_index-1)*N_Timesteps+1;
        if eclipse(current_timestep)==0
            
            %% Measurements
            y_real = real_model.msrFun(x,msrCookieFinal(mag_field_eci(:,current_timestep),...
                sun_pos_eci(:,current_timestep),eclipse(current_timestep),[0;0;0]));
            y_noise = y_real + sqrt(R)*randn(size(y_real));
            [gyro_noise,real_bias] = gyro_noise_func(real_bias,dt,sigma_u,sigma_v);
            

            y_noise(4:6) = y_real(4:6) + gyro_noise;
%             sign=randi([0 1]); 
%             if sign==0 
%                 sign=-1; 
%             end 
            y_noise(7:9) = css_noise(sun_pos_eci(:,current_timestep),x(1:4),xsat_eci(:,current_timestep),albedo(:,current_timestep),lambda);

            if eclipse((cycle_index-1)*N_Timesteps+1)~=0
                y_noise(7:9)=zeros(3,1);
            else
                y_noise(7:9)=y_noise(7:9)/norm(y_noise(7:9)); 
            end
            y_noise(1:3)=y_noise(1:3)/norm(y_noise(1:3)); 

            %% Wahba
            [q_wahba,~]=wahba(y_noise(7:9),y_noise(1:3),sun_pos_eci(:,current_timestep),mag_field_eci(:,current_timestep)*10^(-9));
            mekf.global_state(1:4) = q_wahba';


            %% Bias timer
            bias_init_counter = bias_init_counter + 1;
            % We add zeros for every timestep before the initialization is finished
            if (bias_init_counter < bias_wahba_loops + 1)        
                quat_pos(:,bias_init_counter) = q_wahba;    
                  %Set measurements to 0 until bias has initialized

                for i=1:10
                    %% Current SGP4 matrices values
                    Nodem = nodem(1,current_timestep);
                    Inclm = inclm(1,current_timestep);
                    Argpm = argpm(1,current_timestep);
                    Mm = mm(1,current_timestep);
                    Sun_pos_orbit = sun_pos_orbit(:,current_timestep);
                    Mag_field_orbit = mag_field_orbit(:,current_timestep)*10^(-9);
                    
                    %% Propagate the system
                    current_timestep = (cycle_index-1)*N_Timesteps+i;
                
                    q_ob = quat_EB2OB(x(1:4),Nodem,Inclm,Argpm,Mm);
                    [T_dist, ~,~,~,~,~,~] = disturbances_pd(q_ob, Sun_pos_orbit,Mag_field_orbit, disturbancesEnabled);
                
                    torq = T_dist;

                    x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,rw_ang_momentum,[0;0;0])); 
                    x_real(:,current_timestep)=x;
                    t = t + dt;
                    
                    %% Matrices update
                    x_hat_data(:,current_timestep) =  zeros(7,1);
                    q_ob_data(:,current_timestep) = quat_EB2OB(x(1:4), nodem(1,current_timestep),...
                        inclm(1,current_timestep),argpm(1,current_timestep),mm(1,current_timestep) );
                    bias_data(:,current_timestep) = real_bias;
                    gyro_noise_data(:,current_timestep) = gyro_noise;
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
            mekf.global_state(5:7) = initial_bias_estimate; %Initialize angular velocity equal to gyroscope measurement 

        end
    end 
       %% Main continuous loop
        for cycle_index = cycle_index:number_of_cycles
  
             for timestep_index = 1:3
                 
                current_timestep = (cycle_index-1)*N_Timesteps+timestep_index + 1;
                
                %% Current SGP4 matrices values
                Mag_field_eci = mag_field_eci(:,current_timestep);
                Sun_pos_eci = sun_pos_eci(:,current_timestep);
                Eclipse = eclipse(current_timestep);
                Xsat_eci = xsat_eci(:,current_timestep);
                Albedo_inaccurate = albedo_inaccurate(:,current_timestep);
                Albedo = albedo(:,current_timestep);
                Nodem = nodem(1,current_timestep);
                Inclm = inclm(1,current_timestep);
                Argpm = argpm(1,current_timestep);
                Mm = mm(1,current_timestep);
                Sun_pos_orbit = sun_pos_orbit(:,current_timestep);
                Mag_field_orbit = mag_field_orbit(:,current_timestep)*10^(-9);

                %% Q covariance update
                
                Q_selection(Eclipse,Param.Q,Param.R_hat,mekf,Q_eclipse_load);
                
                %% Sensor Measurements
                y_real = real_model.msrFun(x,msrCookieFinal(Mag_field_eci,Sun_pos_eci,Eclipse,[0;0;0]));
                
                y_noise = y_real + sqrt(R)*randn(size(y_real));
                [gyro_noise,real_bias] = gyro_noise_func(real_bias,dt,sigma_u,sigma_v);

                y_noise(4:6) = y_real(4:6) + gyro_noise;
             
                y_noise(7:9) = css_noise(Sun_pos_eci,x(1:4),Xsat_eci,Albedo,lambda);
        
                if eclipse(current_timestep)~=0
                    y_noise(7:9)=zeros(3,1);
                else
                    y_noise(7:9)=y_noise(7:9)/norm(y_noise(7:9)); 
                end
                y_noise(1:3)=y_noise(1:3)/norm(y_noise(1:3)); 
        
                
                %% MEKF correct
                gyro = y_noise(4:6);
                mekf.correct(y_noise, msrCookieFinalExtended(Mag_field_eci,Sun_pos_eci,Eclipse,gyro,Xsat_eci,Albedo_inaccurate,lambda));

                x_hat = mekf.global_state;
                x_hat(1:4) = x_hat(1:4) / norm(x_hat(1:4));
                
                %% Propagate the system
                              
                [~, ~, ~, AngVel_rw_radps(3,1), acceleration_rw(2,1), T_rw_total] = rw_model(0, AngVel_rw_radps(3,1)); % RW model
                q_ob = quat_EB2OB(x(1:4), Nodem,Inclm,Argpm,Mm );
                [T_dist, ~,~,ad,r,sp,g] = disturbances_pd(q_ob,Sun_pos_orbit, Mag_field_orbit, disturbancesEnabled);
                
                torq = T_rw_total + T_dist;
                x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,rw_ang_momentum,[0;0;0])); 
                  
                
                %% MEKF predict
                % Predict the states at next time step, k+1. This updates the State and
                % StateCovariance properties of the filter to contain x[k+1|k] and
                % P[k+1|k]. These will be utilized by the filter at the next time step.

                gyro = y_noise(4:6);
                mekf.predict(stateTransCookieFinalNominal(torq,rw_ang_momentum,gyro),dt);
                
                residual(timestep_index,:) = y_noise(1:6) - mekf.H_k*mekf.local_error_state;
                %% Matrices update
                x_hat_data(:,current_timestep) = x_hat;
                tau_ad(:,current_timestep) = ad;
                tau_rm(:,current_timestep) = r;
                tau_sp(:,current_timestep) = sp;
                tau_g(:,current_timestep) = g;
                tau_dist(:,current_timestep) = T_dist;                
                x_real(:,current_timestep)=x;
                AngVel_rw_rpm(3,1) = AngVel_rw_radps(3,1)*30/pi; 
                rw_ang_vel_rpm(1,current_timestep) = AngVel_rw_rpm(3,1);  
                Bbody_data(:,current_timestep) = y_real(1:3)*norm(mag_field_orbit(:,current_timestep)*10^(-9));
                bias_data(:,current_timestep) = real_bias;
                gyro_noise_data(:,current_timestep) = gyro_noise;
                q_ob_data(:,current_timestep) = q_ob;
                
                %% Check if the time for Detumbling has come
                
                if current_timestep > 1
                    [trigger_flag, trigger_flag_raw, threshold_times, threshold_exceptions] = ...
                        trigger_N2D(x_real(5:7, current_timestep), x_real(5:7, current_timestep-1), threshold_times, threshold_exceptions);

                    bdot_activation_matrix(1, current_timestep) = trigger_flag;
                    bdot_activation_matrix(2, current_timestep) = trigger_flag_raw;
                end
                
             end 

             for timestep_index=4:10   
        
                 current_timestep = (cycle_index-1)*N_Timesteps+timestep_index +1 ;
                 
                 %% Current SGP4 matrices values 

                 Mag_field_eci = mag_field_eci(:,current_timestep);
                 Sun_pos_eci = sun_pos_eci(:,current_timestep);
                 Eclipse = eclipse(current_timestep);
                 Nodem = nodem(1,current_timestep);
                 Inclm = inclm(1,current_timestep);
                 Argpm = argpm(1,current_timestep);
                 Mm = mm(1,current_timestep);
                 Sun_pos_orbit = sun_pos_orbit(:,current_timestep);
                 Mag_field_orbit = mag_field_orbit(:,current_timestep)*10^(-9);
                 
                 %% Q covariance update
                 
                 Q_selection(Eclipse,Param.Q,Param.R_hat,mekf,Q_eclipse_load);
                
                %% Sensor Measurements
                y_real = real_model.msrFun(x,msrCookieFinal(Mag_field_eci,Sun_pos_eci,Eclipse,[0;0;0]));
                
                [gyro_noise,real_bias] = gyro_noise_func(real_bias,dt,sigma_u,sigma_v);
                y_noise(4:6) = y_real(4:6) + gyro_noise; 

                x_hat = mekf.global_state;
                x_hat(1:4) = x_hat(1:4) / norm(x_hat(1:4));
                
                

                %% PD function
                    
                q_ob_hat = quat_EB2OB(x_hat(1:4),Nodem,Inclm,Argpm,Mm);
                acceleration_rw(1,1) = acceleration_rw(2,1);
                acceleration_rw(2,1) = acceleration_rw(3,1);
                AngVel_rw_radps(1,1) = AngVel_rw_radps(2,1);
                AngVel_rw_radps(2,1) = AngVel_rw_radps(3,1);
                AngVel_rw_rpm(1,1) = AngVel_rw_rpm(2,1);
                AngVel_rw_rpm(2,1) = AngVel_rw_rpm(3,1);

                [torq, T_rw, T_magnetic_effective, ~, ~, ~, AngVel_rw_rpm_next, AngVel_rw_radps_next,...
                        acceleration_rw_cur, rw_ang_momentum, init_AngVel_dz, init_accel_dz, ~, ~, ~, ...
                            timeflag_dz,M] = ...
                                PD(Eclipse,Kp_gain, Kd_gain, q_desired ,q_ob_hat, Const.w_o_io, y_noise(4:6)-mekf.global_state(5:7) , y_noise(1:3)*norm(Mag_field_orbit), ...
                                    Const.mtq_max, Const.lim_dz, AngVel_rw_radps(2,1), AngVel_rw_rpm(2,1), ...
                                        acceleration_rw(1,1), init_AngVel_dz, init_accel_dz, timeflag_dz,Const.rw_max_torque,...
                                            y_real(1:3)*norm(Mag_field_orbit), cycle_index, Const.known_rm);
               
                
                %% Propagate the system
                
                q_ob = quat_EB2OB(x(1:4),Nodem,Inclm,Argpm,Mm);
                [T_dist, ~,~,ad,r,sp,g] = disturbances_pd(q_ob,Sun_pos_orbit, Mag_field_orbit, disturbancesEnabled);
                
                torq = torq + T_dist;
                
                x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,rw_ang_momentum,[0;0;0])); 
                

                %% MEKF predict
                gyro = y_noise(4:6);
                mekf.predict(stateTransCookieFinalNominal(torq,rw_ang_momentum,gyro),dt);
                
                %% Matrices update
                tau_ad(:,current_timestep) = ad;
                tau_rm(:,current_timestep) = r;
                tau_sp(:,current_timestep) = sp;
                tau_g(:,current_timestep) = g;
                tau_dist(:,current_timestep) = T_dist;
                x_real(:,current_timestep)=x;
                tau_rw(1, current_timestep) = T_rw(3);
                tau_mtq(:, current_timestep) = T_magnetic_effective;
                rw_ang_vel_rpm(1,current_timestep) = AngVel_rw_rpm(3,1); 
                acceleration_rw(2,1) = acceleration_rw_cur;
                AngVel_rw_rpm(3,1) = AngVel_rw_rpm_next;
                AngVel_rw_radps(3,1) = AngVel_rw_radps_next;
                M_data(:,current_timestep) = M;
                bias_data(:,current_timestep) = real_bias;
                gyro_noise_data(:,current_timestep) = gyro_noise;
                x_hat_data(:,current_timestep) = x_hat;
                Bbody_data(:,current_timestep) = y_real(1:3)*norm(mag_field_orbit(:,current_timestep)*10^(-9));
                q_ob_data(:,current_timestep) = q_ob;
                
                %% Check if the time for Detumbling has come
                
                if current_timestep > 1
                    [trigger_flag, trigger_flag_raw, threshold_times, threshold_exceptions] = ...
                        trigger_N2D(x_real(5:7, current_timestep), x_real(5:7, current_timestep-1), threshold_times, threshold_exceptions);

                    bdot_activation_matrix(1, current_timestep) = trigger_flag;
                    bdot_activation_matrix(2, current_timestep) = trigger_flag_raw;
                end
                
             end
        end
    
end