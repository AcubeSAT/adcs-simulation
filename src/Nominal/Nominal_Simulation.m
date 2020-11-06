close all;
clc;

clear variables;

set_matlab_utils_path();

%% Initialize Parameters Script
Const=constants();
setParamsFinal_Nominal();
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
t = 0;
n_steps = length(Time)*dt; % Number of time steps
n_steps = uint16(n_steps);
timeflag_dz = 0;
init_AngVel_dz = 0;
init_accel_dz = 0;
rw_ang_momentum=0;
rw_ang_vel_rpm = zeros(1, length(Time));    % RW angular velocity in rpm
rw_accel = zeros(1, length(Time));          % RW acceleration
tau_mtq = zeros(3, length(Time));           % Torques produced by MTQs
tau_rw = zeros(1, length(Time));            % Torques produced by the RW
lambda=1;

%% Initialize matrices
AngVel_rw_radps = zeros(3, 1); %1 = old, 2 = cur, 3 = next
AngVel_rw_rpm = zeros(3, 1);
acceleration_rw = zeros(3, 1);

%% Next we initialize the bias estimation by solving Wahba's problem n times. 

bias_init_counter = 0;  % how many Wahba's have we solved for bias init
bias_wahba_loops = 2;   % Total times to be solved
quat_pos = zeros(4,bias_wahba_loops); % Wahba results are stored here
real_bias=init_bias;

for l=1:n_steps
    
    % Initial state when we get sun sensor measurement
    if eclipse((l-1)/dt+1)==0
        
        %% Measurements
        y_real = real_model.msrFun(x_real(:,(l-1)/dt+1),msrCookieFinal(mag_field_eci(:,(l-1)/dt+1),...
            sun_pos_eci(:,(l-1)/dt+1),eclipse((l-1)/dt+1),[0;0;0]));
        y_noise = y_real + sqrt(R)*randn(size(y_real));
        [gyro_noise,real_bias] = gyro_noise_func(real_bias,dt,sigma_u,sigma_v);
        bias_data = [bias_data real_bias];
        
    %     y_noise=y_real;
        y_noise(4:6) = y_real(4:6) + gyro_noise;
        sign=randi([0 1]); 
        if sign==0 
            sign=-1; 
        end 
        y_noise(7:9) = y_real(7:9) + sign*0.01*y_real(7:9).*poissrnd(lambda,3,1); 
        if eclipse((l-1)/dt+1)~=0
            y_noise(7:9)=zeros(3,1);
        else
            y_noise(7:9)=y_noise(7:9)/norm(y_noise(7:9)); 
        end
        y_noise(1:3)=y_noise(1:3)/norm(y_noise(1:3)); 

    %     y_noise=y_real;

        %% Wahba
        [q_wahba,~]=wahba(y_noise(7:9),y_noise(1:3),sun_pos_eci(:,(l-1)/dt+1),mag_field_eci(:,(l-1)/dt+1));
        ekf.theta(1:4) = q_wahba';


        %% Bias timer
        bias_init_counter = bias_init_counter + 1;
        % We add zeros for every timestep before the initialization is finished
        if (bias_init_counter < bias_wahba_loops + 1)        
            quat_pos(:,bias_init_counter) = q_wahba;    
            x_hat_data = [x_hat_data zeros(7,7)];  %Set measurements to 0 until bias has initialized
    %         y_hat_data = [y_hat_data zeros(n_msr,10)]; 
    %         P_data = [P_data zeros(10*10,10)];

            for i=1:1/dt
                %% Propagate the system
                torq(1,1) = -(-(Const.I(3, 3) - Const.I(2, 2)) * x(6) * x(7))*(1 + 1e-2*randn())-kd*x(5);
                torq(2,1) = -(-(Const.I(1, 1) - Const.I(3, 3)) * x(5) * x(7))*(1 + 1e-2*randn())-kd*x(6);
                torq(3,1) = -(-(Const.I(2, 2) - Const.I(1, 1)) * x(5) * x(6))*(1 + 1e-2*randn())-kd*x(7);

                x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,rw_ang_momentum,[0;0;0])); 
                x_real(:,(l-1)/dt+i+1)=x;
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

     q_ob_data(:,(k-1)/dt+1) = quat_EB2OB(x(1:4), nodem(1,(k-1)/dt+1),inclm(1,(k-1)/dt+1),argpm(1,(k-1)/dt+1),mm(1,(k-1)/dt+1) );   
     for c=1:3
        y_real = real_model.msrFun(x_real(:,(k-1)/dt+c),msrCookieFinal(mag_field_eci(:,(k-1)/dt+c),...
            sun_pos_eci(:,(k-1)/dt+c),eclipse((k-1)/dt+c),[0;0;0]));
        y_noise = y_real + sqrt(R)*randn(size(y_real));
        [gyro_noise,real_bias] = gyro_noise_func(real_bias,dt,sigma_u,sigma_v);
        bias_data = [bias_data real_bias];
        y_noise(4:6) = y_real(4:6) + gyro_noise;
        sign=randi([0 1]); 
        if sign==0 
            sign=-1; 
        end 
        y_noise(7:9) = y_real(7:9) + sign*0.01*y_real(7:9).*poissrnd(lambda,3,1); 
        if eclipse((k-1)/dt+c)~=0
            y_noise(7:9)=zeros(3,1);
        else
            y_noise(7:9)=y_noise(7:9)/norm(y_noise(7:9)); 
        end
        y_noise(1:3)=y_noise(1:3)/norm(y_noise(1:3)); 
        
        gyro = y_noise(4:6);
        ekf.correct(y_noise, msrCookieFinal(mag_field_eci(:,(k-1)/dt+c),...
            sun_pos_eci(:,(k-1)/dt+c),eclipse((k-1)/dt+c),gyro));

        x_hat = ekf.theta;
        x_hat(1:4) = x_hat(1:4) / norm(x_hat(1:4));
        x_hat_data = [x_hat_data x_hat];

        %% Propagate the system
        torq=zeros(3,1);
        [T_dist,~] = disturbances_pd(q_ob_data(:,(k-1)/dt+c), sun_pos_orbit(:,(k-1)/dt+c), mag_field_orbit(:,(k-1)/dt+c)*10^(-9), disturbancesEnabled);
        torq = torq + T_dist;
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,rw_ang_momentum,[0;0;0])); 
        x_real(:,(k-1)/dt+c+1)=x;
        q_ob_data(:,(k-1)/dt+c+1) = quat_EB2OB(x(1:4), nodem(1,(k-1)/dt+c),inclm(1,(k-1)/dt+c),argpm(1,(k-1)/dt+c),mm(1,(k-1)/dt+c) );
        rw_ang_vel_rpm(1,(k-1)/dt+c+1) = AngVel_rw_rpm(3,1); 
    %     if norm(q_prev - q_ob_data(:,(k-1)/dt+c+1))^2 > norm(q_prev + q_ob_data(:,(k-1)/dt+c+1))^2
    %         q_ob_data(:,(k-1)/dt+c+1) = -q_ob_data(:,(k-1)/dt+c+1);
    %     end
%          if q_ob_data(1,(k-1)/dt+c+1)<0
%             q_ob_data(:,(k-1)/dt+c+1) = -q_ob_data(:,(k-1)/dt+c+1);
%          end  
    %    q_prev = q_ob_data(:,(k-1)/dt+c+1);
        % Predict the states at next time step, k+1. This updates the State and
        % StateCovariance properties of the filter to contain x[k+1|k] and
        % P[k+1|k]. These will be utilized by the filter at the next time step.
        
        gyro = y_noise(4:6);
        ekf.predict(stateTransCookieFinalNominal(torq,rw_ang_momentum,gyro));
     end 
    %  
    %  q_ob_hat_prev =quat_EB2OB(x_hat(1:4), nodem(1,(k-1)/dt+c),inclm(1,(k-1)/dt+c),argpm(1,(k-1)/dt+c),mm(1,(k-1)/dt+c) );
    % 
    % if q_ob_hat_prev(1) < 0
    %    q_ob_hat_prev = -q_ob_hat_prev;
    % end
     for c=4:10   
         
        y_real = real_model.msrFun(x_real(:,(k-1)/dt+c),msrCookieFinal(mag_field_eci(:,(k-1)/dt+c),...
            sun_pos_eci(:,(k-1)/dt+c),eclipse((k-1)/dt+c),[0;0;0]));
        [gyro_noise,real_bias] = gyro_noise_func(real_bias,dt,sigma_u,sigma_v);
        bias_data = [bias_data real_bias];
        y_noise(4:6) = y_real(4:6) + gyro_noise; 
        
        x_hat = ekf.theta;
        x_hat(1:4) = x_hat(1:4) / norm(x_hat(1:4));
        x_hat_data = [x_hat_data x_hat];
        q_ob_hat = quat_EB2OB(x_hat(1:4), nodem(1,(k-1)/dt+c-1),inclm(1,(k-1)/dt+c-1),argpm(1,(k-1)/dt+c-1),mm(1,(k-1)/dt+c-1) );

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
                    timeflag_dz] = ...
                        PD(q_desired ,q_ob_hat, ekf.theta(5:7) , y_noise(1:3)*norm(mag_field_orbit(:,(k-1)/dt+c)*10^(-9)) , eclipse((k-1)/dt+c), ...
                            Const.mtq_max, Const.lim_dz, AngVel_rw_radps(2,1), AngVel_rw_rpm(2,1), ...
                                acceleration_rw(1,1), init_AngVel_dz, init_accel_dz, timeflag_dz,Const.rw_max_torque);
        tau_rw(1, (k-1)/dt+c+1) = T_rw(3);
        tau_mtq(:, (k-1)/dt+c+1) = T_magnetic_effective;
        rw_ang_vel_rpm(1,(k-1)/dt+c+1) = AngVel_rw_rpm(3,1); 
        acceleration_rw(2,1) = acceleration_rw_cur;
        AngVel_rw_rpm(3,1) = AngVel_rw_rpm_next;
        AngVel_rw_radps(3,1) = AngVel_rw_radps_next;
        
        %%                    
        [T_dist,~] = disturbances_pd(q_ob_data(:,(k-1)/dt+c), sun_pos_orbit(:,(k-1)/dt+c), mag_field_orbit(:,(k-1)/dt+c)*10^(-9), disturbancesEnabled);
        torq = torq + T_dist;
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,rw_ang_momentum,[0;0;0])); 
        x_real(:,(k-1)/dt+c+1)=x;
        q_ob_data(:,(k-1)/dt+c+1) = quat_EB2OB(x(1:4), nodem(1,(k-1)/dt+c),inclm(1,(k-1)/dt+c),argpm(1,(k-1)/dt+c),mm(1,(k-1)/dt+c) );



    %     if norm(q_prev - q_ob_data(:,(k-1)/dt+c+1))^2 > norm(q_prev + q_ob_data(:,(k-1)/dt+c+1))^2
    %         q_ob_data(:,(k-1)/dt+c+1) = -q_ob_data(:,(k-1)/dt+c+1);
    %     end
%         if q_ob_data(1,(k-1)/dt+c+1)<0
%             q_ob_data(:,(k-1)/dt+c+1) = -q_ob_data(:,(k-1)/dt+c+1);
%         end  
    %     q_prev = q_ob_data(:,(k-1)/dt+c+1);
        gyro = y_noise(4:6);
        ekf.predict(stateTransCookieFinalNominal(torq,rw_ang_momentum,gyro));
     end
        end
        break;
    else
        for i=1:1/dt
            %% Propagate the system
            torq=zeros(3,1);

            x = real_model.stateTransFun(x, stateTransCookieFinal(torq)); 
            x_real(:,(l-1)/dt+i+1)=x;
            t = t + dt;

        end
        x_hat_data = [x_hat_data zeros(10,10)];
    end

end
%% Calculation of performance error
x_real_euler_perf = quat2eul(q_ob_data(1:4,1:length(q_ob_data)-1)'); 
x_real_euler_perf = rad2deg(x_real_euler_perf');

instant_error_perform = x_real_euler_perf';

figure();
for i=1:3
    subplot(3,1,i);
    hold on;
    plot(Time(21:length(instant_error_perform)), instant_error_perform(21:length(instant_error_perform), i), 'LineWidth',1.5, 'Color','blue');
    if (i==1), title('Absolute Performance Errors', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis'); end
    if (i==2), ylabel('Y-axis'); end
    if (i==3), ylabel('Z-axis'); end
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

%% Calulation of knowledge error
x_hat_euler_know = zeros(length(x_hat_data), 6);
instant_error_know = zeros(length(x_hat_data), 6);

x_real_euler_know = quat2eul(x_real(1:4,1:length(x_hat_data))');
x_real_euler_know = rad2deg(x_real_euler_know');
x_hat_euler_know(:, 1:3) = quat2eul(x_hat_data(1:4,:)');
x_hat_euler_know(:, 1:3) = (rad2deg(x_hat_euler_know(:, 1:3)'))';

instant_error_know(:, 1:3) = x_hat_euler_know(:, 1:3) - x_real_euler_know';
instant_error_know(:, 4:6) = x_hat_data(5:7, 1:length(x_hat_data))' - x_real(5:7, 1:length(x_hat_data))';


for k=1:length(instant_error_know)
    if instant_error_know(k, 1) > 180
        instant_error_know(k, 1) = instant_error_know(k, 1) - 180;
    elseif instant_error_know(k, 1) < -180
        instant_error_know(k, 1) = instant_error_know(k, 1) + 180;
    end
end

figure();
for i=1:6
    subplot(6,1,i);
    hold on;
    plot(Time(21:length(instant_error_know)), instant_error_know(21:length(instant_error_know), i), 'LineWidth',1.5, 'Color','blue');
    ylabel(['$\tilde{x}_' num2str(i) '$'], 'interpreter','latex', 'fontsize',14);
    if (i==1), title('Absolute Knowledge Errors', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis'); end
    if (i==2), ylabel('Y-axis'); end
    if (i==3), ylabel('Z-axis'); end
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

n_dim = size(x_real,1);
% figure('Position',[500 0 1420 1080]);
figure();
for i=1:n_dim
    subplot(n_dim,1,i);
    hold on;
    plot(Time,x_real(i,1:length(Time)), 'LineWidth',2.0, 'Color','blue');
    plot(Time(1:length(x_hat_data(i,:))),x_hat_data(i,:), 'LineWidth',2.0, 'Color','magenta');
    if (i==1),legend({['$x_' num2str(i) '$'],['$\hat{x}_' num2str(i) '$']}, 'interpreter','latex', 'fontsize',15);end
    ylabel(['$x_' num2str(i) '$'], 'interpreter','latex', 'fontsize',17);
    if (i==1), title('EKF estimation results', 'interpreter','latex', 'fontsize',17);end
    xlim([3 n_steps]);
    hold off;
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
    hold on;
    plot(Time,q_ob_data(i,1:length(Time)), 'LineWidth',2.0, 'Color','blue');
    ylabel(['$q_{ob' num2str(i) '}$'], 'interpreter','latex', 'fontsize',17);
    xlim([3 n_steps]);
    hold off;
end

figure()
plot(1:length(eclipse),eclipse, 'LineWidth',2.0, 'Color','blue');
xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
ylabel(['Eclipse'], 'interpreter','latex', 'fontsize',14);
if (i==1), title('Umbral, Penumbral or no Eclipse', 'interpreter','latex', 'fontsize',17);end

x_err_data = x_real(:,1:length(x_hat_data))-x_hat_data(1:7,:);
figure();
for i=1:n_dim
    subplot(n_dim,1,i);
    plot(Time(1:length(x_err_data(i,:))),x_err_data(i,:), 'LineWidth',2.0, 'Color','blue');  % Error for the first state
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    ylabel(['$\tilde{x}_' num2str(i) '$'], 'interpreter','latex', 'fontsize',14);
    if (i==1), legend({'State estimate','$\pm \sigma$'}, 'interpreter','latex', 'fontsize',11);end
    if (i==1), title('State estimation errors', 'interpreter','latex', 'fontsize',11); end
    xlim([3 n_steps]);
    hold off;
end

eul_diff =zeros(length(Time),3);
euler_hat=quat2eul(q_ob_data(1:4,:)')';
eul_diff=(euler_hat)*180/pi;
figure();
for i=1:3
    subplot(3,1,i);
    hold on;
    plot(Time,eul_diff(i,1:length(Time)), 'LineWidth',2.0, 'Color','blue');
    ylabel(['$euler_{' num2str(i) '}$'], 'interpreter','latex', 'fontsize',17);
    xlim([3 n_steps]);
    hold off;
end
%%  Plotting the produced Torques
 
      figure()
      subplot(3,1,1)
      plot(1:length(Time),tau_mtq(1,1:length(Time)))
      title('Magnetic Torques-x')
      ylabel('Torque [Nm]')
      xlabel('Time [s]')
      grid on;
      subplot(3,1,2)
      plot(1:length(Time),tau_mtq(2,1:length(Time)))
      title('Magnetic Torques-y')
      ylabel('Torque [Nm]')
      xlabel('Time [s]')
      grid on;
      subplot(3,1,3)
      plot(1:length(Time),tau_mtq(3,1:length(Time)))
      title('Magnetic Torques-z')
      ylabel('Torque [Nm]')
      xlabel('Time [s]')
      grid on;
    
      figure()
      plot(1:length(Time),tau_rw(1, 1:length(Time)))
      title('Reaction Wheel Torque-z')
      ylabel('Torque [Nm]')
      xlabel('Time [s]')
      grid on;

 figure() 
 plot(Time, rw_ang_vel_rpm(1,1:length(Time)),'LineWidth',1.5, 'Color','blue'); 
 title('Angular velocity of RW'); 
 ylabel('Angular Velocity [rpm]'); 
 xlabel('Time [s]'); 
 grid on; 