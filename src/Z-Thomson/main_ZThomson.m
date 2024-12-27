% ========================================================================
%   Main function for Z-Thomson Mode simulation.
% ========================================================================

close all;
clc;
clear variables;

%% Initialize Parameters Script
Const = constants();
Param = setParams_ZThomson(Const.I);

dt_model = Param.dt_model;
orbits = Param.orbits;
x0 = Param.x0;
x0_hat = Param.x0_hat;
real_model = Param.real_model;
model = Param.model;

albedo = Param.albedo;
albedo_inaccurate = Param.albedo_inaccurate;
setDisturbances = Param.setDisturbances;

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
total_limit = Param.total_limit;
cycle_duration = Param.cycle_duration;
Total_simulation_time = Param.tst;
R_old = Param.R_old;
satrec = Param.satrec;
Kd = Param.Kd;
Ks = Param.Ks;
w_ref = Param.w_ref;

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
n_params = length(x0_hat);                                                          % Init number of parameters
mekf = MEKF(n_params, number_of_measurements, @model.stateTransFun, @model.msrFun); % Init EKF
mekf.global_state = x0_hat;                                                         % Init state estimation
mekf.P = P0;                                                                        % Init Covariance matrix
mekf.setProcessNoiseCov(Q);                                                         % Q variance matrix
mekf.setMeasureNoiseCov(R_hat);                                                     % R variance matrix
mekf.setFadingMemoryCoeff(1.00);                                                    % This parameter defined the "memory" of the filter, 1 meaning regular
mekf.setPartDerivStep(0.001);                                                       % Arithmetical jacobian step
if (use_analytic_jacob)                                                             % If analytical jacobian is used, define the functions
    mekf.setStateTransFunJacob(@model.stateTransFunJacob);
    mekf.setMsrFunJacob(@model.msrFunJacob);
end

%% Initialize simulation parameters

Time = 0:cycle_duration:Total_simulation_time;
x_real = zeros(7, length(Time)); % Real state
q_ob_data = zeros(4,length(Time));
x = x0(1:7);
x_real(:,1) = x0(1:7);
t = 0;
B_body = zeros(3, length(Time));                    % Earth magnetic field satellite in body frame
w_b_ob = zeros(3, length(Time));                    % Angular velocity from orbit to body frame expressed in body frame
w_b_ib = zeros(3,length(Time));
w_o_ob = zeros(3,length(Time));
Bdot_body = zeros(3, length(Time));                 % Bdot estimation in body frane
Bdot_body_z = zeros(1, length(Time));
Bdot_body_new = zeros(3, length(Time)); 
w_b_ob_magn = zeros(1, length(Time));               % Angular velocity magnitude
Mag = zeros(3, length(Time));                       % Commanded magnetic dipole moment
w_b_ob_Bdot = zeros(3, length(Time));               % Angular velocity estimation using Bdot
theta_deg_arr= zeros(3,length(Time));
z_body_in_orbit = zeros(3,length(Time));
theta_deg_arr_x = zeros(1,length(Time));
theta_deg_arr_y = zeros(1,length(Time));
theta_deg_arr_z = zeros(1,length(Time));
q_orbit_body_data = zeros(4,length(Time));
w=zeros(3,length(Time));


tau_mtq = zeros(3, length(Time));           % Torques produced by MTQs
M_data = zeros(3, length(Time));
tau_dist = zeros(3, length(Time));
lambda=1;
tau_ad = zeros(3, length(Time));
tau_rm = zeros(3, length(Time));
tau_g = zeros(3, length(Time));
tau_sp = zeros(3, length(Time));
plotter_step = 1;
reps = length(Time);

x_hat_data = zeros(7,length(Time));
bias_data = zeros(3,length(Time));
gyro_noise_data = zeros(3,length(Time));
Bbody_data = zeros(3,length(Time));

bias_wahba_loops = 2;
N_Timesteps = 10;
for cycle_index = 1:bias_wahba_loops
    current_timestep = (cycle_index - 1) * N_Timesteps + 1;
    if eclipse(current_timestep) == 0

        %% Measurements

        y_real = real_model.msrFun(x, msrCookieFinal(mag_field_eci(:, current_timestep), ...
            sun_pos_eci(:, current_timestep), eclipse(current_timestep), [0; 0; 0]));
        y_noise = y_real + sqrt(R) * randn(size(y_real));
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



for current_cycle = 1:length(Time) %Main loop

    %% Angular velocity and Magnetic field expressed in body frame
    q_orbit_eci = dcm2quat(Orbit2ECI_DCM(nodem(1, current_cycle), inclm(1, current_cycle), argpm(1, current_cycle)+mm(1, current_cycle)));
    q_eci_body = x(1:4);
    q_orbit_body = quatProd(q_orbit_eci, q_eci_body);
    q_orbit_body_data(:,current_cycle) = q_orbit_body;
    R_OB = quat2dcm(q_orbit_body'); % Calculating the transformation matrix from orbit to body frame
    B_body(:, current_cycle) = R_OB * mag_field_orbit(:, current_cycle) * 10^(-9);
    B_body(:, current_cycle) = B_body(:, current_cycle) + sqrt(R_old) * randn(size(B_body(:, current_cycle))); % Adding white noise to magnetometer measurements
    R_BO = R_OB';
    w_b_io = R_OB(:, 3) * Const.w_o;
    w_b_ob(:, current_cycle) = x(5:7) - w_b_io; % Calculating angular rate of satellite relative to ECI frame
    w_b_ob_magn(current_cycle) = norm(w_b_ob(:, current_cycle));
    w_b_ib(:,current_cycle)=  w_b_ob(:, current_cycle) + w_b_io;
    w_o_ob(:,current_cycle)= R_BO*w_b_ob(:,current_cycle);
    w(:,current_cycle)=rotate_vector([0 ;1 ;1; 1],w_b_ib(:,current_cycle));
    %% First timestep - Orbit propagation
    [T_disturbances, ~] = disturbances_bdot(R_BO, sun_pos_orbit(:, current_cycle), B_body(:, current_cycle), setDisturbances); % Calculation of external torques
    for j = 1:(cycle_duration / dt_model * 0.1)
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(T_disturbances, 0, 0));
    end
    [~, ~, ~, ~, ~, ~, ~, ~, ~, mag_field_orbit2, ~, ~, ~, ~, ~, argpm2, nodem2, inclm2, mm2, ~, ~] = orbit_sgp4_offset(satrec, cycle_duration, cycle_duration, current_cycle-1+.1);

    %% Angular velocity and Magnetic field expressed in body frame
    q_orbit_eci = dcm2quat(Orbit2ECI_DCM(nodem2, inclm2, argpm2+mm2));
    q_eci_body = x(1:4);
    q_orbit_body = quatProd(q_orbit_eci, q_eci_body);
    R_OB = quat2dcm(q_orbit_body'); % Calculating the transformation matrix from orbit to body frame
    B_body_2 = R_OB * mag_field_orbit2 * 10^(-9);
    B_body_2 = B_body_2 + sqrt(R_old) * randn(size(B_body(:, current_cycle))); % Second measurment from magnetometer
    R_EB=quat2dcm(q_eci_body');
    %% Bdot calculation
    Bx=acos(B_body(1,current_cycle)/norm(B_body(:,current_cycle)));
    Bx2=acos(B_body_2(1)/norm(B_body_2));

    By=acos(B_body(2,current_cycle)/norm(B_body(:,current_cycle)));
    By2=acos(B_body_2(2)/norm(B_body_2));
    
    Bz=acos(B_body(3,current_cycle)/norm(B_body(:,current_cycle)));
    Bz2=acos(B_body_2(3)/norm(B_body_2));

    Bdot_body(:, current_cycle) = (B_body_2 - B_body(:, current_cycle)) / 0.1;
    Bdot_body_new(:, current_cycle) = ([Bx2;By2;Bz2] - [Bx;By;Bz]) / 0.1;
    Bdot_body_z(current_cycle) = (Bz2 - Bz) / 0.1;

    %Bdot_body(1:2, current_cycle) = (B_body_2(1:2) - B_body(1:2, current_cycle)) / 0.1;
    %Bdot_body(3, current_cycle) = (Bz2 - Bz) / 0.1;

    %% Torque Calculation
    %% Z-Thomson using 3 mtqs
    %M(1,1) = -Kd*(w_b_ob_Bdot(3,current_cycle)-0.1)/norm(w_b_ob_Bdot(3,current_cycle)-0.1)*sign(Bdot_body(2,current_cycle));
    %M(2,1)= -Kd*(w_b_ob_Bdot(3,current_cycle)-0.1)/norm(w_b_ob_Bdot(3,current_cycle)-0.1)*sign(Bdot_body(1,current_cycle));
    %M(3,1) = -Kd * Bdot_body(3, current_cycle) / norm(B_body(3, current_cycle));
    %M = mtq_scaling(M, Const.mtq_max); % MTQ scaling
    %Mag(:, current_cycle) = M;
    %T_magnetic = cross(M, B_body(:, current_cycle));
    
    
    %% Z-Thomson using 2 mtqs
    % M(1,1) = -Kd*(w_b_ob_Bdot(3,current_cycle)-0.1)/norm(w_b_ob_Bdot(3,current_cycle)-0.1)*sign(Bdot_body(2,current_cycle));
    % M(2,1)= 0;
    % M(3,1) = -Kd * Bdot_body(3, current_cycle) / norm(B_body(3, current_cycle));
    % M = mtq_scaling(M, Const.mtq_max); % MTQ scaling
    % Mag(:, current_cycle) = M;
    % T_magnetic = cross(M, B_body(:, current_cycle));
 
   
      
      % %% Z Thomson New attempt
      % 
      % M(1,1) = -Ks * (w_b_oo(3,current_cycle) - 0.05) * sign(B_body(2,current_cycle));
      % %M(1,1) = 0;
      % %M(2,1) = Ks * (w_b_oo(3,current_cycle) - 0.05) * sign(B_body(1,current_cycle));
      % M(2,1) = 0;
      % %M(3,1) = -Kd * Bdot_body(3,current_cycle) / norm(B_body(:,current_cycle)); % the old way
      % M(3,1) = Kd * Bdot_body_z(current_cycle);  % the new new way 
      % 
      % M = mtq_scaling(M, Const.mtq_max); % MTQ scaling
      % Mag(:, current_cycle) = M;
      % T_magnetic = cross(M, B_body(:, current_cycle));
        
       %% Z Thomson based on lappas' paper

       % if (abs(w_b_ib) < 0.12)
       %     w_ref = w_ref2;
       % else
       %     w_ref = w_ref1;
       % end

         % M(3,1)= Kd * Bdot_body_z(current_cycle);
         % if (abs(B_body(2,current_cycle))>abs(B_body(1,current_cycle)))
         % M(1,1)= Ks*(w_b_ib(3,current_cycle)-w_ref)*sign(B_body(2,current_cycle));
         % M(2,1)= 0;
         % elseif (abs(B_body(2,current_cycle))<abs(B_body(1,current_cycle)))
         %  M(1,1)= 0;
         %  M(2,1)=-Ks*(w_b_ib(3,current_cycle)-w_ref)*sign(B_body(1,current_cycle));
         % end  

         M(3,1)= Kd * Bdot_body_z(current_cycle);
         if (abs(B_body(2,current_cycle))>abs(B_body(1,current_cycle)))
            M(1,1)= -Ks*(abs(w_b_ib(3,current_cycle))-w_ref)*sign(B_body(2,current_cycle));
            M(2,1)= 0;
         elseif (abs(B_body(2,current_cycle))<abs(B_body(1,current_cycle)))
            M(1,1)= 0;
            M(2,1)=Ks*(abs(w_b_ib(3,current_cycle))-w_ref)*sign(B_body(1,current_cycle));
         end  


         M = mtq_scaling(M, Const.mtq_max); % MTQ scaling
         Mag(:, current_cycle) = M;
         T_magnetic = cross(M, B_body(:, current_cycle));

        %% Angle between Z-Axis of b.f. and Z-Axis of o.f.
        x_ob=R_OB(:,1); 
        y_ob=R_OB(:,2);
        z_ob=R_OB(:,3); % Extract the third column of R_OB(z-axis of orbit in body frame)
        z_body=[0; 0; 1]; % z-axis of the body frame
        cos_theta_x=dot(x_ob,z_body)/(norm(x_ob)*norm(z_body)); % Calculate cosine of the angle
        cos_theta_y=dot(y_ob,z_body)/(norm(y_ob)*norm(z_body));
        cos_theta_z=dot(z_ob,z_body)/(norm(z_ob)*norm(z_body));
                     
        theta_x=acos(cos_theta_x); % compute angle in radians
        theta_y=acos(cos_theta_y);
        theta_z=acos(cos_theta_z);
        theta_deg_x=rad2deg(theta_x); %Convert to degrees
        theta_deg_y=rad2deg(theta_y) ;
        theta_deg_z=rad2deg(theta_z) ;
        theta_deg_arr_x(current_cycle) = theta_deg_x;
        theta_deg_arr_y(current_cycle) = theta_deg_y;
        theta_deg_arr_z(current_cycle) = theta_deg_z;

    %% Calculate V, I, P of MTQ's
    [V_mtq, I_mtq, P_thermal_mtq] = mtq_model(M);

    %% Estimation of angular velocity using Bdot
    w_b_ob_Bdot(:, current_cycle) = skew(B_body(:, current_cycle)) * (-Bdot_body_new(:, current_cycle) ...
        -B_body(:, current_cycle)) / (B_body(:, current_cycle)' * B_body(:, current_cycle));



    %%  Disturbances for 2 timesteps to acquire magnetometer measurements
   [T_disturbances, ~] = disturbances_bdot(R_BO, sun_pos_orbit(:, current_cycle), B_body(:, current_cycle), setDisturbances);
   % T_disturbances=[0;0;0];
   torq = T_magnetic + T_disturbances;
    for j = 1:(cycle_duration / dt_model * 0.2)
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(T_disturbances, 0, 0));
    end

    %% Apply calculated torque for last 7 timesteps
    for j = 1:(cycle_duration / dt_model * 0.7)
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq, 0, 0));
    end
    t = t + cycle_duration;

end







%%  Plotting the Angular Velocities
figure()
subplot(3, 1, 1)
plot(1:plotter_step:length(Time), w_b_ob(1, 1:plotter_step:end))
title('Angular Velocities w_b_ob', 'interpreter', 'latex', 'fontsize', 17);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel(['$\omega_', num2str(1), '$', '[rad/sec]'], 'interpreter', 'latex', 'fontsize', 14);
grid on
subplot(3, 1, 2)
plot(1:plotter_step:length(Time), w_b_ob(2, 1:plotter_step:end))
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel(['$\omega_', num2str(2), '$', '[rad/sec]'], 'interpreter', 'latex', 'fontsize', 14);
grid on
subplot(3, 1, 3)
plot(1:plotter_step:length(Time), w_b_ob(3, 1:plotter_step:end))
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel(['$\omega_', num2str(3), '$', '[rad/sec]'], 'interpreter', 'latex', 'fontsize', 14);
grid on

%%  Plotting the Angular Velocities
figure()
subplot(3, 1, 1)
plot(1:plotter_step:length(Time), w_b_ib(1, 1:plotter_step:end))
title('Angular Velocities w_b_ib', 'interpreter', 'latex', 'fontsize', 17);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel(['$\omega_', num2str(1), '$', '[rad/sec]'], 'interpreter', 'latex', 'fontsize', 14);
grid on
subplot(3, 1, 2)
plot(1:plotter_step:length(Time), w_b_ib(2, 1:plotter_step:end))
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel(['$\omega_', num2str(2), '$', '[rad/sec]'], 'interpreter', 'latex', 'fontsize', 14);
grid on
subplot(3, 1, 3)
plot(1:plotter_step:length(Time), w_b_ib(3, 1:plotter_step:end))
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel(['$\omega_', num2str(3), '$', '[rad/sec]'], 'interpreter', 'latex', 'fontsize', 14);
grid on

%%  Plotting the Angular Velocities
figure()
subplot(3, 1, 1)
plot(1:plotter_step:length(Time), w_o_ob(1, 1:plotter_step:end))
title('Angular Velocities w_o_ob', 'interpreter', 'latex', 'fontsize', 17);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel(['$\omega_', num2str(1), '$', '[rad/sec]'], 'interpreter', 'latex', 'fontsize', 14);
grid on
subplot(3, 1, 2)
plot(1:plotter_step:length(Time), w_o_ob(2, 1:plotter_step:end))
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel(['$\omega_', num2str(2), '$', '[rad/sec]'], 'interpreter', 'latex', 'fontsize', 14);
grid on
subplot(3, 1, 3)
plot(1:plotter_step:length(Time), w_o_ob(3, 1:plotter_step:end))
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel(['$\omega_', num2str(3), '$', '[rad/sec]'], 'interpreter', 'latex', 'fontsize', 14);
grid on


% Plotting the angular velocity using Bdot
figure()
subplot(3, 1, 1)
plot(1:plotter_step:length(Time), w_b_ob_Bdot(1, 1:plotter_step:end))
title('Angular Velocity using Bdot')
xlabel('Time [s]');
ylabel('Angular Velocity-x [rad/sec]');
grid on;
subplot(3, 1, 2)
plot(1:plotter_step:length(Time), w_b_ob_Bdot(2, 1:plotter_step:end))
xlabel('Time [s]');
ylabel('Angular Velocity-y [rad/sec]');
grid on;
subplot(3, 1, 3)
plot(1:plotter_step:length(Time), w_b_ob_Bdot(3, 1:plotter_step:end))
xlabel('Time [s]');
ylabel('Angular Velocity-z [rad/sec]');
grid on;







 % Bdot
% figure()
% subplot(3, 1, 1)
% plot(1:plotter_step:length(Time), Bdot_body(1, 1:plotter_step:end))
% title('Bdot-x');
% xlabel('Time[s]');
% ylabel('Bdot [T/sec)]');
% subplot(3, 1, 2)
% plot(1:plotter_step:length(Time), Bdot_body(2, 1:plotter_step:end))
% title('Bdot-y');
% xlabel('Time[s]');
% ylabel('Bdot [T/sec)]');
% subplot(3, 1, 3)
% plot(1:plotter_step:length(Time), Bdot_body(3, 1:plotter_step:end))
% title('Bdot-z');
% xlabel('Time[s]');
% ylabel('Bdot [T/sec)]');

%  Plotting the Angular Velocity Magnitude
figure()
plot(1:plotter_step:length(Time), w_b_ob_magn(1:plotter_step:end))
title('Angular Velocity Magnitude');
xlabel('Time [s]');
ylabel('|Velocity| [rad/sec]');

%  Plotting the produced dipole moment
figure()
subplot(3, 1, 1)
plot(1:plotter_step:length(Time), Mag(1, 1:plotter_step:end))
title('Magnetic Dipole-x');
xlabel('Time [s]');
ylabel({'Magnetic'; 'Dipole [Am^2]'});
subplot(3, 1, 2)
plot(1:plotter_step:length(Time), Mag(2, 1:plotter_step:end))
title('Magnetic Dipole-y');
xlabel('Time [s]');
ylabel({'Magnetic'; 'Dipole [Am^2]'});
subplot(3, 1, 3)
plot(1:plotter_step:length(Time), Mag(3, 1:plotter_step:end))
title('Magnetic Dipole-z');
xlabel('Time [s]');
ylabel({'Magnetic'; 'Dipole [Am^2]'});



% Plotting eclipse
% figure()
% plot(1:length(eclipse), eclipse, 'LineWidth', 2.0, 'Color', 'blue');
% title('Eclipse', 'interpreter', 'latex', 'fontsize', 17)
% xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
% ylabel('Eclipse', 'interpreter', 'latex', 'fontsize', 14);
% grid on;
% if (i == 1), title('Umbral, Penumbral or no Eclipse', 'interpreter', 'latex', 'fontsize', 17); end


% Plotting the Angle between Z-Axis of b.f. and Z-Axis of o.f.
% figure()
% plot(1:plotter_step:length(Time), theta_deg_arr(1,1:plotter_step:end));
% title('Angle between X-Axis of b.f. and Z-Axis of o.f.', 'interpreter', 'latex', 'fontsize', 17);
% xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
% ylabel('Angle [degrees]', 'interpreter', 'latex', 'fontsize', 14);
% grid on
% 
% figure()
% plot(1:plotter_step:length(Time), theta_deg_arr(2,1:plotter_step:end));
% title('Angle between Y-Axis of b.f. and Z-Axis of o.f.', 'interpreter', 'latex', 'fontsize', 17);
% xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
% ylabel('Angle [degrees]', 'interpreter', 'latex', 'fontsize', 14);
% grid on
% 
% figure()
% plot(1:plotter_step:length(Time), theta_deg_arr(3,1:plotter_step:end));
% title('Angle between Z-Axis of b.f. and Z-Axis of o.f.', 'interpreter', 'latex', 'fontsize', 17);
% xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
% ylabel('Angle [degrees]', 'interpreter', 'latex', 'fontsize', 14);
% grid on

% Plotting the Angle between Z-Axis of b.f. and Z-Axis of o.f.
figure()
plot(1:plotter_step:length(Time), theta_deg_arr_x(1:plotter_step:end));
title('Angle between Z-Axis of b.f. and X-Axis of o.f.', 'interpreter', 'latex', 'fontsize', 17);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Angle [degrees]', 'interpreter', 'latex', 'fontsize', 14);
grid on

figure()
plot(1:plotter_step:length(Time), theta_deg_arr_y(1:plotter_step:end));
title('Angle between Z-Axis of b.f. and Y-Axis of o.f.', 'interpreter', 'latex', 'fontsize', 17);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Angle [degrees]', 'interpreter', 'latex', 'fontsize', 14);
grid on


figure()
plot(1:plotter_step:length(Time), theta_deg_arr_z(1:plotter_step:end));
title('Angle between Z-Axis of b.f. and Z-Axis of o.f.', 'interpreter', 'latex', 'fontsize', 17);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Angle [degrees]', 'interpreter', 'latex', 'fontsize', 14);
grid on

figure()
plot(B_body(3,:))
title('Magnetic field in Z axis', 'interpreter', 'latex', 'fontsize', 17);

figure()
plot(B_body(2,:))
title('Magnetic field in Y axis', 'interpreter', 'latex', 'fontsize', 17);

figure()
plot(B_body(1,:))
title('Magnetic field in X axis', 'interpreter', 'latex', 'fontsize', 17);

% figure()
% grid on;
% hold on
% for i = 1:length(Time)
%     % Clear current plot
%     cla;
% 
%     % Plot vector as an arrow (from origin [0, 0, 0] to [x(t), y(t), z(t)])
%     quiver3(0, 0, 0, z_body_in_orbit(1,i), z_body_in_orbit(2,i), z_body_in_orbit(3,i), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
% 
%     % Plot the trajectory of the tip of the vector
%     plot3(z_body_in_orbit(1,1:i), z_body_in_orbit(2,1:i),z_body_in_orbit(3,1:i), 'b--'); % trajectory of vector time
% 
%     % Plot the vector as an arrow (from origin [0, 0, 0] to [x(t), y(t), z(t)])
%     %quiver3(0, 0, 0, z_body_in_orbit(1,i), z_body_in_orbit(2,i), z_body_in_orbit(3,i), 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
% 
%     % Plot the trajectory of the tip of the vector (the path it traces)
%     %plot3(z_body_in_orbit(1,1:i), z_body_in_orbit(2,1:i),z_body_in_orbit(3,1:i), 'b--', 'LineWidth', 1.5); % Blue dashed line
% 
%     % Set axis limits and labels
%     xlim([-1.2, 1.2]);
%     ylim([-1.2, 1.2]);
%     zlim([-1.2, 1.2]);
% 
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
% 
%     % Pause for animation effect
%     pause(0.1);
% end
% 
% hold off;

% Parameters for the 3D rectangle
width = 1;   % Width along the X-axis
height = 1;  % Height along the Y-axis
depth = 3;   % Depth along the Z-axis

% Define the rectangular prism's vertices
vertices = [ -width/2, -height/2, -depth/2;
              width/2, -height/2, -depth/2;
              width/2,  height/2, -depth/2;
             -width/2,  height/2, -depth/2;
             -width/2, -height/2,  depth/2;
              width/2, -height/2,  depth/2;
              width/2,  height/2,  depth/2;
             -width/2,  height/2,  depth/2];

% Define the faces for patch plotting
faces = [1 2 3 4;
         5 6 7 8;
         1 2 6 5;
         2 3 7 6;
         3 4 8 7;
         4 1 5 8];

% Initialize the figure and axis
figure;
axis equal;
grid on;
xlim([-5 5]); ylim([-5 5]); zlim([-5 5]);
xlabel('X Orbit Frame'); ylabel('Y Orbit Frame'); zlabel('Z Orbit Frame');
hold on;

% Set the view for 3D visualization
view(3);  % This ensures the plot uses a 3D view


% Create patch object for the rectangular prism
rect_prism = patch('Vertices', vertices, 'Faces', faces, ...
                   'FaceColor', 'cyan', 'EdgeColor', 'black');

% Initialize text handle for the current time and frame
time_text = text(0, 5, 5, '', 'FontSize', 12, 'Color', 'black');  % Position the text above the plot area

% Animation loop
for t = 1:length(Time)

    % Rotation matrix
    q_orbit_body = q_orbit_body_data(:,t);
    R_OB = quat2dcm(q_orbit_body');

    % Apply rotation to each vertex
    rotated_vertices = (R_OB' * vertices')';  % Transform vertices using R

    % Update the rectangular prism's vertices
    set(rect_prism, 'Vertices', rotated_vertices);

    % Update the time and frame label
    set(time_text, 'String', sprintf('Time: %d', t));

    % Update the plot
    drawnow;

    pause(0.005);
end