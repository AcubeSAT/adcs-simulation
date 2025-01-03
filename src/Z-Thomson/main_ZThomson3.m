% ========================================================================
%   Main function for Z Thomson Mode simulation.
% ========================================================================

close all;
clear all;
clc;

%% Initialize Parameters Script

Const = constants();
Param = setParams_ZThomson3(Const.I);
satrec=Param.satrec;
dt = Param.dt;
orbits = Param.orbits;
tf = Param.tf;
q_desired = Param.q_desired;
x0 = Param.x0;
x0_hat = Param.x0_hat;
real_model = Param.real_model;
model = Param.model;
number_of_measurements = Param.number_of_measurements;
albedo = Param.albedo;
albedo_inaccurate = Param.albedo_inaccurate;
disturbancesEnabled = Param.disturbancesEnabled;
xsat_eci = Param.xsat_eci;
eclipse = Param.eclipse;
mag_field_eci = Param.mag_field_eci;
mag_field_orbit = Param.mag_field_orbit;
sun_pos_eci = Param.sun_pos_eci;
sun_pos_orbit = Param.sun_pos_orbit;
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
use_analytic_jacob = Param.use_analytic_jacob;
total_limit = Param.total_limit;
exceptions_limit= Param.exceptions_limit;
N_Timesteps= Param.N_Timesteps;


%% Initialize Global Parameters

global R_coils; R_coils = Const.R_coils;
global N_coils; N_coils = Const.N_coils;
global A_coils; A_coils = Const.A_coils;
global mt; mt = Const.mt;


%% Construct the MEKF

Q_struct = load('idealQ.mat', 'IDEALQ');
Q_eclipse_load = Q_struct.IDEALQ;
n_params = length(x0_hat);                                                              % Init number of parameters
mekf = MEKF(n_params, number_of_measurements, @model.stateTransFun, @model.msrFun);     % Init EKF
mekf.global_state = x0_hat;                                                             % Init state estimation
mekf.P = P0;                                                                            % Init Covariance matrix
mekf.setProcessNoiseCov(Q);                                                             % Q variance matrix
mekf.setMeasureNoiseCov(R_hat);                                                         % R variance matrix
mekf.setFadingMemoryCoeff(1.00);                                                        % This parameter defined the "memory" of the filter, 1 meaning regular
mekf.setPartDerivStep(0.001);                                                           % Arithmetical jacobian step
if (use_analytic_jacob)                                                                 % If analytical jacobian is used, define the functions
    mekf.setStateTransFunJacob(@model.stateTransFunJacob);
    mekf.setMsrFunJacob(@model.msrFunJacob);
end

%% Initialize simulation parameters

Time = 0:dt:tf;
x_real = zeros(7, length(Time)); % Real state
q_ob_data = zeros(4, length(Time));
x = x0(1:7);
x_real(:, 1) = x0(1:7);
t = 0;
number_of_cycles = floor(length(Time)/N_Timesteps); % Number of cycles

rw_ang_momentum = 0;

tau_mtq = zeros(3, length(Time));                   % Torques produced by MTQs
M_data = zeros(3, length(Time));
tau_dist = zeros(3, length(Time));
lambda = 1;
tau_ad = zeros(3, length(Time));
tau_rm = zeros(3, length(Time));
tau_g = zeros(3, length(Time));
tau_sp = zeros(3, length(Time));
plotter_step = 1;
reps = length(Time);

x_hat_data = zeros(7, length(Time));
bias_data = zeros(3, length(Time));
gyro_noise_data = zeros(3, length(Time));
Bbody_data = zeros(3, length(Time));
q_sb_data = zeros(4, length(Time));

estimated_velocity = zeros(3, length(Time));
theta_deg_arr_x= zeros(3, length(Time));
theta_deg_arr_y= zeros(3, length(Time));
theta_deg_arr_z= zeros(3, length(Time));
y_noise = zeros(9,1);

%% Next we initialize the bias estimation by solving Wahba's problem n times.

bias_init_counter = 0; % how many Wahba's have we solved for bias init
bias_wahba_loops = 2; % Total times to be solved
quat_pos = zeros(4, bias_wahba_loops); % Wahba results are stored here
real_bias = init_bias;


for cycle_index = 1:bias_wahba_loops
    current_timestep = (cycle_index - 1) * N_Timesteps + 1;
    if eclipse(current_timestep) == 0
        
        Mag_field_orbit = mag_field_orbit(:,current_timestep)*10^(-9);

        %% Measurements
        y_real = real_model.msrFun(x, msrCookieFinal(mag_field_eci(:, current_timestep), ...
            sun_pos_eci(:, current_timestep), eclipse(current_timestep), [0; 0; 0]));
        y_noise(1:3) = y_real(1:3)*norm(Mag_field_orbit) + 1e-9*diag([15,15,15])*randn(3,1);
        [gyro_noise, real_bias] = gyro_noise_func(real_bias, dt, sigma_u, sigma_v);


        y_noise(4:6) = y_real(4:6) + gyro_noise;
        % sign=randi([0 1]);
        % if sign==0
        %     sign=-1;
        % end
        y_noise(7:9) = css_noise(sun_pos_eci(:, current_timestep), x(1:4), xsat_eci(:, current_timestep), albedo(:, current_timestep), lambda);

        if eclipse((cycle_index - 1)*N_Timesteps+1) ~= 0
            y_noise(7:9) = zeros(3, 1);
        else
            y_noise(7:9) = y_noise(7:9) / norm(y_noise(7:9));
        end
        y_noise(1:3) = y_noise(1:3) / norm(y_noise(1:3));

        %% Wahba
        [q_wahba, ~] = wahba(y_noise(7:9), y_noise(1:3), sun_pos_eci(:, current_timestep), mag_field_eci(:, current_timestep)*10^(-9));
        mekf.global_state(1:4) = q_wahba';

        %% Bias timer
        
        bias_init_counter = bias_init_counter + 1;
        quat_pos(:, bias_init_counter) = q_wahba;
        % We add zeros for every timestep before the initialization is finished
        if (bias_init_counter < bias_wahba_loops + 1)
            
            %Set measurements to 0 until bias has initialized

            for i = 1:10

                %% Current SGP4 matrices values
                    Nodem = nodem(1,current_timestep);
                    Inclm = inclm(1,current_timestep);
                    Argpm = argpm(1,current_timestep);
                    Mm = mm(1,current_timestep);
                    Sun_pos_orbit = sun_pos_orbit(:,current_timestep);
                    Mag_field_orbit = mag_field_orbit(:,current_timestep)*10^(-9);

                %% Propagate the system
                current_timestep = (cycle_index - 1) * N_Timesteps + i;

                q_ob = quat_EB2OB(x(1:4), Nodem, Inclm, Argpm, Mm);
                [T_dist, ~, ~, ~, ~, ~, ~] = disturbances_pd(q_ob, Sun_pos_orbit, Mag_field_orbit, disturbancesEnabled);

                torq = T_dist;
                
                x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq, rw_ang_momentum, [0;0;0]));
                x_real(:, current_timestep) = x;
                t = t + dt;

                %% Matrices update
                    x_hat_data(:,current_timestep) =  zeros(7,1);
                    q_ob_data(:,current_timestep) = quat_EB2OB(x(1:4), nodem(1,current_timestep),...
                        inclm(1,current_timestep),argpm(1,current_timestep),mm(1,current_timestep) );
                    bias_data(:,current_timestep) = real_bias;
                    gyro_noise_data(:,current_timestep) = gyro_noise;
                     R_OB = quat2dcm(q_ob');
                 % Angle between Z-Axis of b.f. and Z-Axis of o.f.
                     x_ob=R_OB(:,1) ;
                    y_ob=R_OB(:,2);
                    z_ob=R_OB(:,3); % Extract the third column of R_OB(z-axis of orbit in body frame)
                    z_body=[0; 0; 1]; % z-axis of the body frame
                    cos_theta_x=dot(x_ob,z_body)/(norm(x_ob)*norm(z_body)); % Calculate cosine of the angle
                    cos_theta_y=dot(y_ob,z_body)/(norm(y_ob)*norm(z_body));
                    cos_theta_z=dot(z_ob,z_body)/(norm(z_ob)*norm(z_body));
            
                    theta_x=acos(cos_theta_x); % compute angle in radians
                    theta_y=acos(cos_theta_y);
                    theta_z=acos(cos_theta_z);
                    theta_x=rad2deg(theta_x); %Convert to degrees
                    theta_y=rad2deg(theta_y) ;
                    theta_z=rad2deg(theta_z) ;
                    theta_deg_arr_x(:, current_timestep)=theta_x; %Convert to degrees
                    theta_deg_arr_y(:, current_timestep)=theta_y ;
                    theta_deg_arr_z(:, current_timestep)=theta_z ;

            end
%             continue
        end

        %% Bias calculation
        dqdt = zeros(4, bias_wahba_loops-1);
        for i = 2:bias_wahba_loops
            dqdt(:, i-1) = quat_pos(:, i) - quat_pos(:, i-1);
        end

        estimated_rate = zeros(3, bias_wahba_loops-1);

        % The angular rate is calculated using: angular_rate = 2 * dq/dt * q^-1
        for i = 1:bias_wahba_loops - 1
            temp = 2 * quatProd(quatconj(quat_pos(:, i)'), dqdt(:, i));
            estimated_rate(:, i) = temp(2:4);
        end

        real_rate = y_noise(4:6);

        mean_omega = [0; 0; 0]; % Time averaging for noise reduction

        for i = 1:3
            mean_omega(i) = mean(estimated_rate(i, :));
        end

        initial_bias_estimate = real_rate - mean_omega;
        mekf.global_state(5:7) = initial_bias_estimate; %Initialize angular velocity equal to gyroscope measurement

    end
end

%% Main continuous loop

for cycle_index = cycle_index:number_of_cycles
    for timestep_index = 1:3

        current_timestep = (cycle_index - 1) * N_Timesteps + timestep_index + 1;

        %% Current SGP4 matrices values
        Mag_field_eci = mag_field_eci(:, current_timestep);
        Sun_pos_eci = sun_pos_eci(:, current_timestep);
        Eclipse = eclipse(current_timestep);
        Xsat_eci = xsat_eci(:, current_timestep);
        Albedo_inaccurate = albedo_inaccurate(:, current_timestep);
        Albedo = albedo(:, current_timestep);
        Nodem = nodem(1, current_timestep);
        Inclm = inclm(1, current_timestep);
        Argpm = argpm(1, current_timestep);
        Mm = mm(1, current_timestep);
        Sun_pos_orbit = sun_pos_orbit(:, current_timestep);
        Mag_field_orbit = mag_field_orbit(:, current_timestep) * 10^(-9);

        %% Q covariance update
%         Q_selection(Eclipse, Param.Q, Param.R_hat, mekf, Q_eclipse_load);

        %% Sensor Measurements
        y_real = real_model.msrFun(x, msrCookieFinal(Mag_field_eci, Sun_pos_eci, Eclipse, [0; 0; 0]));

        y_noise(1:3) = y_real(1:3)*norm(Mag_field_orbit) + 1e-9*diag([15,15,15])*randn(3,1);
        [gyro_noise, real_bias] = gyro_noise_func(real_bias, dt, sigma_u, sigma_v);

        y_noise(4:6) = y_real(4:6) + gyro_noise;

        y_noise(7:9) = css_noise(Sun_pos_eci, x(1:4), Xsat_eci, Albedo, lambda);

        if eclipse(current_timestep)~=0
            y_noise(7:9) = zeros(3, 1);
        else
            y_noise(7:9) = y_noise(7:9) / norm(y_noise(7:9));
        end
        y_noise(1:3) = y_noise(1:3) / norm(y_noise(1:3));

        %% MEKF correct
        gyro = y_noise(4:6);
        mekf.correct(y_noise, msrCookieFinalExtended(Mag_field_eci, Sun_pos_eci, Eclipse, gyro, Xsat_eci, Albedo_inaccurate, lambda));

        x_hat = mekf.global_state;
        x_hat(1:4) = x_hat(1:4) / norm(x_hat(1:4));

        %% Propagate the system

        q_ob = quat_EB2OB(x(1:4), Nodem, Inclm, Argpm, Mm);
        [T_dist, ~, ~, ad, r, sp, g] = disturbances_pd(q_ob, Sun_pos_orbit, Mag_field_orbit, disturbancesEnabled);

        torq =  T_dist;
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq, 0, [0;0;0]));

        %% MEKF predict
        % Predict the states at next time step, k+1. This updates the State and
        % StateCovariance properties of the filter to contain x[k+1|k] and
        % P[k+1|k]. These will be utilized by the filter at the next time step.

        gyro = y_noise(4:6);
        mekf.predict(stateTransCookieFinalNominal(torq, 0, gyro), dt);

        %% Matrices update
        x_hat_data(:, current_timestep) = x_hat;
        tau_ad(:, current_timestep) = ad;
        tau_rm(:, current_timestep) = r;
        tau_sp(:, current_timestep) = sp;
        tau_g(:, current_timestep) = g;
        tau_dist(:, current_timestep) = T_dist;
        x_real(:, current_timestep) = x;

        Bbody_data(:, current_timestep) = y_real(1:3) * norm(mag_field_orbit(:, current_timestep)*10^(-9));


        bias_data(:, current_timestep) = real_bias;
        gyro_noise_data(:, current_timestep) = gyro_noise;
        q_ob_data(:, current_timestep) = q_ob;
        q_sb_data(:, current_timestep) = q_sun_body(Sun_pos_eci, x(1:4),Const.sun_desired)';
        R_OB = quat2dcm(q_ob');
                % Angle between Z-Axis of b.f. and Z-Axis of o.f.
        x_ob=R_OB(:,1) ;
        y_ob=R_OB(:,2);
        z_ob=R_OB(:,3); % Extract the third column of R_OB(z-axis of orbit in body frame)
        z_body=[0; 0; 1]; % z-axis of the body frame
        cos_theta_x=dot(x_ob,z_body)/(norm(x_ob)*norm(z_body)); % Calculate cosine of the angle
        cos_theta_y=dot(y_ob,z_body)/(norm(y_ob)*norm(z_body));
        cos_theta_z=dot(z_ob,z_body)/(norm(z_ob)*norm(z_body));

        theta_x=acos(cos_theta_x); % compute angle in radians
        theta_y=acos(cos_theta_y);
        theta_z=acos(cos_theta_z);
        theta_x=rad2deg(theta_x); %Convert to degrees
        theta_y=rad2deg(theta_y) ;
        theta_z=rad2deg(theta_z) ;
        theta_deg_arr_x(:, current_timestep)=theta_x; %Convert to degrees
        theta_deg_arr_y(:, current_timestep)=theta_y ;
        theta_deg_arr_z(:, current_timestep)=theta_z ;
        sun_orbit_normalized = (Sun_pos_orbit / norm(Sun_pos_orbit));
        Const.sun_desired = Const.sun_desired / norm(Const.sun_desired);

        estimated_velocity(:, current_timestep) = gyro - x_hat(5:7);

        if timestep_index == 2
            B_body_thomson = y_noise(1:3)*norm(Mag_field_orbit);
        end
        
    end

    for timestep_index = 4:10

        current_timestep = (cycle_index - 1) * N_Timesteps + timestep_index + 1;

        %% Current SGP4 matrices values
        Mag_field_eci = mag_field_eci(:, current_timestep);
        Sun_pos_eci = sun_pos_eci(:, current_timestep);
        Eclipse = eclipse(current_timestep);
        Nodem = nodem(1, current_timestep);
        Inclm = inclm(1, current_timestep);
        Argpm = argpm(1, current_timestep);
        Mm = mm(1, current_timestep);
        Sun_pos_orbit = sun_pos_orbit(:, current_timestep);
        Mag_field_orbit = mag_field_orbit(:, current_timestep) * 10^(-9);

        %% Q covariance update
%         Q_selection(Eclipse, Param.Q, Param.R_hat, mekf, Q_eclipse_load);

        %% Sensor Measurements
        y_real = real_model.msrFun(x, msrCookieFinal(Mag_field_eci, Sun_pos_eci, Eclipse, [0; 0; 0]));

        [gyro_noise, real_bias] = gyro_noise_func(real_bias, dt, sigma_u, sigma_v);
        y_noise(4:6) = y_real(4:6) + gyro_noise;

        x_hat = mekf.global_state;
        x_hat(1:4) = x_hat(1:4) / norm(x_hat(1:4));

        %% PD function
        % Choose x_hat for determination, x for ground truth 

        q_ob_hat = quat_EB2OB(x_hat(1:4),Nodem,Inclm,Argpm,Mm);
        % q_ob_hat = quat_EB2OB(x(1:4),Nodem,Inclm,Argpm,Mm);

          q_ob = quat_EB2OB(x(1:4), Nodem, Inclm, Argpm, Mm);
          R_OB=quat2dcm(q_ob');

        % Choose first PD for determination, second PD for ground truth

%           [~, ~, ~, ~, ~, ~, ~, ~, ~, mag_field_orbit2, ~, ~, ~, ~, ~, argpm2, nodem2, inclm2, mm2, ~, ~] = orbit_sgp4_offset(satrec, 1, 1, .1);


         [Bdot_body,torq,V_mtq, I_mtq, P_thermal_mtq,M] = ...
        Control_Thomson(x(5:7), y_noise(1:3)*norm(Mag_field_orbit), B_body_thomson, Const.mtq_max);



        %% Propagate the system
        q_ob = quat_EB2OB(x(1:4), Nodem, Inclm, Argpm, Mm);
        [T_dist, ~, ~, ad, r, sp, g] = disturbances_pd(q_ob, Sun_pos_orbit, Mag_field_orbit, disturbancesEnabled);

        torq = torq + T_dist;
        
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,0,[0;0;0]));

        %% MEKF predict
        gyro = y_noise(4:6);
        mekf.predict(stateTransCookieFinalNominal(torq, 0, gyro), dt);

        %% Matrices update
        tau_ad(:, current_timestep) = ad;
        tau_rm(:, current_timestep) = r;
        tau_sp(:, current_timestep) = sp;
        tau_g(:, current_timestep) = g;
        tau_dist(:, current_timestep) = T_dist;
        x_real(:, current_timestep) = x;





        M_data(:, current_timestep) = M;
        bias_data(:, current_timestep) = real_bias;
        gyro_noise_data(:, current_timestep) = gyro_noise;
        x_hat_data(:, current_timestep) = x_hat;
        Bbody_data(:, current_timestep) = y_real(1:3) * norm(mag_field_orbit(:, current_timestep)*10^(-9));
        q_ob_data(:, current_timestep) = q_ob;
        q_sb_data(:, current_timestep) = q_sun_body(Sun_pos_eci, x(1:4),Const.sun_desired)';
        R_OB = quat2dcm(q_ob');
                 % Angle between Z-Axis of b.f. and Z-Axis of o.f.
        x_ob=R_OB(:,1) ;
        y_ob=R_OB(:,2);
        z_ob=R_OB(:,3); % Extract the third column of R_OB(z-axis of orbit in body frame)
        z_body=[0; 0; 1]; % z-axis of the body frame
        cos_theta_x=dot(x_ob,z_body)/(norm(x_ob)*norm(z_body)); % Calculate cosine of the angle
        cos_theta_y=dot(y_ob,z_body)/(norm(y_ob)*norm(z_body));
        cos_theta_z=dot(z_ob,z_body)/(norm(z_ob)*norm(z_body));

        theta_x=acos(cos_theta_x); % compute angle in radians
        theta_y=acos(cos_theta_y);
        theta_z=acos(cos_theta_z);
        theta_x=rad2deg(theta_x); %Convert to degrees
        theta_y=rad2deg(theta_y) ;
        theta_z=rad2deg(theta_z) ;
           theta_deg_arr_x(:, current_timestep)=theta_x; %Convert to degrees
        theta_deg_arr_y(:, current_timestep)=theta_y ;
        theta_deg_arr_z(:, current_timestep)=theta_z ;

        sun_orbit_normalized = (Sun_pos_orbit / norm(Sun_pos_orbit));
        Const.sun_desired = Const.sun_desired / norm(Const.sun_desired);

        estimated_velocity(:, current_timestep) = gyro - x_hat(5:7);


    end
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
        for cycle_index=1:length(instant_error_know)
            if instant_error_know(cycle_index, i) > 180
                instant_error_know(cycle_index, i) = instant_error_know(cycle_index, i) - 360;
            elseif instant_error_know(cycle_index, i) < -180
                instant_error_know(cycle_index, i) = instant_error_know(cycle_index, i) + 360;
            end
        end
    end

    figure();
    for i=1:6
        subplot(6,1,i);
        hold on;
        plot(Time(1:length(instant_error_know)), instant_error_know(1:length(instant_error_know), i), 'LineWidth',1.5, 'Color','blue');
        if (i==1), title('Absolute Knowledge Errors', 'interpreter','latex', 'fontsize',17);end
        if (i==1), ylabel('X-axis [deg]', 'interpreter','latex', 'fontsize',14); end
        if (i==2), ylabel('Y-axis [deg]', 'interpreter','latex', 'fontsize',14); end
        if (i==3), ylabel('Z-axis [deg]', 'interpreter','latex', 'fontsize',14); end
        if (i==4), ylabel('$\omega_1 [rad/sec]$', 'interpreter','latex', 'fontsize',14); end
        if (i==5), ylabel('$\omega_2 [rad/sec]$', 'interpreter','latex', 'fontsize',14); end
        if (i==6), ylabel('$\omega_3 [rad/sec]$', 'interpreter','latex', 'fontsize',14); end
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
        %     xlim([3 number_of_cycles]);
        xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
        hold off;
        grid on;
    end
    
%% Eclipse plot
figure()
plot(1:length(eclipse), eclipse, 'LineWidth', 2.0, 'Color', 'blue');
title('Eclipse', 'interpreter', 'latex', 'fontsize', 17)
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Eclipse', 'interpreter', 'latex', 'fontsize', 14);
grid on;
if (i == 1), title('Umbral, Penumbral or no Eclipse', 'interpreter', 'latex', 'fontsize', 17); end


figure();
for i = 1:3
    subplot(3, 1, i);
    plot(Time(1:end-1), x_real(4+i, 1:length(Time)-1), 'LineWidth', 2.0, 'Color', 'blue');
    xlabel('Time [s]', 'interpreter', 'latex', 'fontsize', 12);
    if (i == 1), ylabel('$\omega_1$ [rad/sec]', 'interpreter', 'latex', 'fontsize', 14); end
    if (i == 2), ylabel('$\omega_2$ [rad/sec]', 'interpreter', 'latex', 'fontsize', 14); end
    if (i == 3), ylabel('$\omega_3$ [rad/sec]', 'interpreter', 'latex', 'fontsize', 14); end
    ylabel(['$\omega_', num2str(i), '$', '[rad/sec]'], 'interpreter', 'latex', 'fontsize', 14);
    if (i == 1), legend('Angular Velocity'); end
    if (i == 1), title('Angular Velocities', 'interpreter', 'latex', 'fontsize', 17); end
    grid on;
end


% Plotting the Angle between Z-Axis of b.f. and Z-Axis of o.f.
figure()
plot(Time(1:end-1), theta_deg_arr_x(1:length(Time)-1));
title('Angle between Z-Axis of b.f. and X-Axis of o.f.', 'interpreter', 'latex', 'fontsize', 17);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Angle [degrees]', 'interpreter', 'latex', 'fontsize', 14);
grid on

figure()
plot(Time(1:end-1), theta_deg_arr_y(1:length(Time)-1));
title('Angle between Z-Axis of b.f. and Y-Axis of o.f.', 'interpreter', 'latex', 'fontsize', 17);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Angle [degrees]', 'interpreter', 'latex', 'fontsize', 14);
grid on


figure()
plot(Time(1:end-1), theta_deg_arr_z(1:length(Time)-1));
title('Angle between Z-Axis of b.f. and Z-Axis of o.f.', 'interpreter', 'latex', 'fontsize', 17);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Angle [degrees]', 'interpreter', 'latex', 'fontsize', 14);
grid on




% % Parameters for the 3D rectangle
% width = 1;   % Width along the X-axis
% height = 1;  % Height along the Y-axis
% depth = 3;   % Depth along the Z-axis
% 
% % Define the rectangular prism's vertices
% vertices = [ -width/2, -height/2, -depth/2;
%               width/2, -height/2, -depth/2;
%               width/2,  height/2, -depth/2;
%              -width/2,  height/2, -depth/2;
%              -width/2, -height/2,  depth/2;
%               width/2, -height/2,  depth/2;
%               width/2,  height/2,  depth/2;
%              -width/2,  height/2,  depth/2];
% 
% % Define the faces for patch plotting
% faces = [1 2 3 4;
%          5 6 7 8;
%          1 2 6 5;
%          2 3 7 6;
%          3 4 8 7;
%          4 1 5 8];
% 
% % Initialize the figure and axis
% figure;
% axis equal;
% grid on;
% xlim([-5 5]); ylim([-5 5]); zlim([-5 5]);
% xlabel('X Orbit Frame'); ylabel('Y Orbit Frame'); zlabel('Z Orbit Frame');
% hold on;
% 
% % Set the view for 3D visualization
% view(3);  % This ensures the plot uses a 3D view
% 
% 
% % Create patch object for the rectangular prism
% rect_prism = patch('Vertices', vertices, 'Faces', faces, ...
%                    'FaceColor', 'cyan', 'EdgeColor', 'black');
% 
% % Initialize text handle for the current time and frame
% time_text = text(0, 5, 5, '', 'FontSize', 12, 'Color', 'black');  % Position the text above the plot area
% 
% % Animation loop
% for t = 1:length(Time)
% 
%     % Rotation matrix
%     q_orbit_body = q_orbit_body_data(:,t);
%     R_OB = quat2dcm(q_orbit_body');
% 
%     % Apply rotation to each vertex
%     rotated_vertices = (R_OB' * vertices')';  % Transform vertices using R
% 
%     % Update the rectangular prism's vertices
%     set(rect_prism, 'Vertices', rotated_vertices);
% 
%     % Update the time and frame label
%     set(time_text, 'String', sprintf('Time: %d', t));
% 
%     % Update the plot
%     drawnow;
% 
%     pause(0.005);
% end