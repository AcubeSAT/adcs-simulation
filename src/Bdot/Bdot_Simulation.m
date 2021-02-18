close all;
clc;

clear variables;

set_matlab_utils_path();

% Initialize Parameters Script
Const=constants();
setParamsFinal_Detumbling();

Time = 0:dt:tf;
x_data = zeros(7,length(Time)); % Real state
x = x0;
t = 0;

% x_data = zeros(7,length(Time));
B_body_ctrl = zeros(3,length(Time));
% torq = zeros(3,length(Time));
w_b_ob = zeros(3,length(Time));
accel_b_ob = zeros(3,length(Time));
Bdot_body = zeros(3,length(Time));
Bdot_body_new = zeros(3, length(Time));
w_b_ob_magn = zeros(1,length(Time));
Mag = zeros(3,length(Time));
tau_mtq = zeros(3, length(Time));
q_error = zeros(4, length(Time));
w_b_ob_Bdot = zeros(3, length(Time)); 
accel_b_ob_Bdot = zeros(3,length(Time));
tau_dist = zeros(3,length(Time));
w_b_ob_Bdot_magn = zeros(1,length(Time));
nominal_activation_matrix = zeros(2, length(Time)); % 1 - Bdot, 2 - nonBdot
threshold_times = 0; 
threshold_exceptions = 0;
%N2D_threshold = 0.087; % rad/sec 
D2N_threshold = 0.035; % rad/sec

%% Disturbances Once & Correct Control Cycle
for i=1:length(Time)


%     [~, ~,~,~, ~,~, ~,~,~,mag_field_orbit1, ~,~,~,~,~,~,~,~,~,~,~] = orbit_sgp4_offset(satrec,dt,dt,i-1);

    q_orbit_eci=dcm2quat(Orbit2ECI_DCM(nodem(1,i),inclm(1,i),argpm(1,i)+mm(1,i)));

    %% w_b_ob/B_body based on Orbit2ECI
    q_eci_body = x(1:4);
    q_orbit_body = quatProd(q_orbit_eci,q_eci_body);

    q_error(:, i) = q_orbit_body;
    
    R_OB = quat2dcm(q_orbit_body'); % Calculating the transformation matrix from orbit to body frame
    B_body_ctrl(:,i)=R_OB*mag_field_orbit(:,i)*10^(-9);

    B_body_ctrl(:,i)=B_body_ctrl(:,i)+sqrt(R)*randn(size(B_body_ctrl(:,i)));
    R_BO = R_OB';

    w_b_io = R_OB(:,3)*Const.w_o; 
    w_b_ob(:,i) = x(5:7) - w_b_io; % Calculating angular rate of satellite relative to ECI frame
    if i > 1
        accel_b_ob(:,i) = (w_b_ob(:,i) - w_b_ob(:,i-1))/1;
    end
    w_b_ob_magn(i) = norm(w_b_ob(:,i));
    [T_disturbances, disturbancesMessage] = disturbances_bdot(R_BO, quat2eul(q_orbit_body'), ... 
      sun_pos_orbit(:,i), B_body_ctrl(:,i), setDisturbances);
for j=1:(dt/dt_model*0.1)
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(T_disturbances,0,0)); 
end
    
    [~, ~,~,~, ~,~, ~,~,~,mag_field_orbit2, ~,~,~,~,~,argpm2,nodem2,inclm2,mm2,~,~] =...
        orbit_sgp4_offset(satrec,dt,dt,i-1+.1);

    q_orbit_eci=dcm2quat(Orbit2ECI_DCM(nodem2,inclm2,argpm2+mm2));

    %% w_b_ob/B_body based on Orbit2ECI
    q_eci_body = x(1:4);
    q_orbit_body = quatProd(q_orbit_eci,q_eci_body);

    R_OB = quat2dcm(q_orbit_body'); % Calculating the transformation matrix from orbit to body frame
    B_body_ctrl2=R_OB*mag_field_orbit2*10^(-9);

    B_body_ctrl2=B_body_ctrl2+sqrt(R)*randn(size(B_body_ctrl(:,i)));
    
    %% Bdot
        Bdot_body(:,i) = (B_body_ctrl2 - B_body_ctrl(:,i))/0.1;
        Bdot_body_new(:,i) = cross(B_body_ctrl(:,i), w_b_ob(:,i));

    %% Torque Calculation
%         rm = [0.05; 0.05; 0.05];
        M = -Kp.*Bdot_body(:,i)./norm(B_body_ctrl(:,i));    % Calculate dipole moment to produce
%         M = M - rm;
        M = mtq_scaling(M, Const.mtq_max);                             % MTQ scaling 
        Mag(:, i) = M;
        T_magnetic = cross(M,B_body_ctrl(:,i));
        
    %% Calculate V, I, P of MTQ's

        %[V_mtq, I_mtq, P_thermal_mtq] = mtq_model(M);    

    %% Calculate the angular velocity using the B-dot 
        %w_b_ob_Bdot(:, i) = -Kp.*(Bdot_body(:,i)'*Bdot_body(:,i))/M; 
        w = skew(B_body_ctrl(:,i))*(-Bdot_body(:, i)-B_body_ctrl(:,i))/(B_body_ctrl(:,i)'*B_body_ctrl(:,i)); 
        %w = skew(B_body_ctrl(:,i))*(-Bdot_body(:, i)+B_body_ctrl(:,i))/(B_body_ctrl(:,i)'*B_body_ctrl(:,i)); 
        %w = skew(B_body_ctrl(:,i))*(-Bdot_body(:, i))/(B_body_ctrl(:,i)'*B_body_ctrl(:,i)); 
 
        w_b_ob_Bdot(:, i) = w; 
        
        if i > 1
            accel_b_ob_Bdot(:,i) = (w_b_ob_Bdot(:,i) - w_b_ob_Bdot(:,i-1))/1;
        end
        w_b_ob_Bdot_magn(i) = norm(w_b_ob_Bdot(:,i)); 
    
    %% Calculate when Nominal is ready to be activated
        
        bdot_activation = 0;
        nonBdot_activation = 0;
        total_limit = 520;
        exceptions_limit = 30;
        
        if i > 1 
            if abs(w_b_ob_Bdot(1,i)) < D2N_threshold ... 
                    && abs(w_b_ob_Bdot(2,i)) < D2N_threshold ... 
                        && abs(w_b_ob_Bdot(3,i)) < D2N_threshold 
                 if threshold_times == 0 
                    threshold_times = 1; 
                 end 
                 if threshold_times >= 1    
                     if abs(w_b_ob_Bdot(1,i-1)) < D2N_threshold ... 
                        && abs(w_b_ob_Bdot(2,i-1)) < D2N_threshold ... 
                            && abs(w_b_ob_Bdot(3,i-1)) < D2N_threshold 
                        threshold_times = threshold_times + 1;
                     elseif (abs(w_b_ob_Bdot(1,i-1)) >= D2N_threshold ... 
                             || abs(w_b_ob_Bdot(2,i-1)) >= D2N_threshold ... 
                                || abs(w_b_ob_Bdot(3,i-1)) >= D2N_threshold) ...
                                    && threshold_exceptions < exceptions_limit
                        threshold_times = threshold_times + 1;
                        threshold_exceptions = threshold_exceptions + 1;
                     elseif (abs(w_b_ob_Bdot(1,i-1)) >= D2N_threshold ... 
                             || abs(w_b_ob_Bdot(2,i-1)) >= D2N_threshold ... 
                                || abs(w_b_ob_Bdot(3,i-1)) >= D2N_threshold) ...
                                    && threshold_exceptions >= exceptions_limit
                        threshold_times = 0;
                        threshold_exceptions = 0; 
                     end 
                 end 
                 if threshold_times >= total_limit
                     bdot_activation = 1;
                 end     
            end 
        end
        
        if i > 1 
            if abs(w_b_ob(1,i)) < D2N_threshold ... 
                    && abs(w_b_ob(2,i)) < D2N_threshold ... 
                        && abs(w_b_ob(3,i)) < D2N_threshold 
                    
                    nonBdot_activation = 1;
                    
            end
        end
        
        nominal_activation_matrix(1,i) = bdot_activation;
        nominal_activation_matrix(2,i) = nonBdot_activation;
        
    %%
    [T_disturbances, disturbancesMessage] = disturbances_bdot(R_BO, quat2eul(q_orbit_body'), ... 
      sun_pos_orbit(:,i), B_body_ctrl(:,i), setDisturbances);
    torq = T_magnetic + T_disturbances; 
    for j=1:(dt/dt_model*0.1)
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(T_disturbances,0,0)); 
    end

    for j=1:(dt/dt_model*0.7)
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,0,0)); 
    end
    for j=1:(dt/dt_model*0.1)
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(T_disturbances,0,0)); 
    end
    t = t + dt;
end

%% Disturbances Once
% for i=1:length(Time)
% %     x_data(:,i) = x;
% 
%     q_orbit_eci=dcm2quat(Orbit2ECI_DCM(nodem(1,i),inclm(1,i),argpm(1,i)+mm(1,i)));
% 
%     %% w_b_ob/B_body based on Orbit2ECI
%     q_eci_body = x(1:4);
%     q_orbit_body = quatProd(q_orbit_eci,q_eci_body);
% 
%     R_OB = quat2dcm(q_orbit_body'); % Calculating the transformation matrix from orbit to body frame
% 
% %     [~, ~,~,~, ~,~, ~,~,~,mag_field_orbit1, ~,~,~,~,~,~,~,~,~,~,~] = orbit_sgp4_offset(satrec,dt,dt,i-1);
% 
% %     [~, ~,~,~, ~,~, ~,~,~,mag_field_orbit2, ~,~,~,~,~,~,~,~,~,~,~] = orbit_sgp4_offset(satrec,dt,dt,i-1+.1);
% 
%     B_body_ctrl(:,i)=R_OB*mag_field_orbit(:,i)*10^(-9);
% 
%     B_body_ctrl(:,i)=B_body_ctrl(:,i)+sqrt(R)*randn(size(B_body_ctrl(:,i)));
%     R_BO = R_OB';
% 
%     w_b_io = R_OB(:,3)*Const.w_o; 
%     w_b_ob(:,i) = x(5:7) - w_b_io; % Calculating angular rate of satellite relative to ECI frame
%     w_b_ob_magn(i) = norm(w_b_ob(:,i));
%     %% Bdot
%     if i==1
%         Bdot_body(:,i)=[0;0;0];
%         M=zeros(3,1);
%         Mag(:,i)= M;
%     elseif i>1
%         Bdot_body(:,i) = (B_body_ctrl(:,i) - B_body_ctrl(:,i-1))/dt;
% 
%         %% Torque Calculation
%         M = -Kp.*Bdot_body(:,i)./norm(B_body_ctrl(:,i));    % Calculate dipole moment to produce
%         M = mtq_scaling(M, mtq_max);                             % MTQ scaling 
%         Mag(:,i)= M;
%     end
%         T_magnetic = cross(M,B_body_ctrl(:,i));
% 
%     [T_disturbances, disturbancesMessage] = disturbances(R_BO, quat2eul(q_orbit_body'), ... 
%       sun_pos_orbit(:,i), B_body_ctrl(:,i), setDisturbances);
%     torq = T_magnetic + T_disturbances; 
%     for j=1:(dt/dt_model*0.3)
%         x = real_model.stateTransFun(x, stateTransCookieFinal(T_disturbances)); 
%     end
% 
%     for j=1:(dt/dt_model*0.7)
%         x = real_model.stateTransFun(x, stateTransCookieFinal(torq)); 
%     end
%     t = t + dt;
% end

    %% Disturbances in each timestep
% for i=1:length(Time)
% %     x_data(:,i) = x;
% 
%     q_orbit_eci=dcm2quat(Orbit2ECI_DCM(nodem(1,(i-1)/dt_model+1),inclm(1,(i-1)/dt_model+1),argpm(1,(i-1)/dt_model+1)+mm(1,(i-1)/dt_model+1)));
% 
%     %% w_b_ob/B_body based on Orbit2ECI
%     q_eci_body = x(1:4);
%     q_orbit_body = quatProd(q_orbit_eci,q_eci_body);
% 
%     R_OB = quat2dcm(q_orbit_body'); % Calculating the transformation matrix from orbit to body frame
% 
%     B_body_ctrl(:,i)=R_OB*mag_field_orbit(:,(i-1)/dt_model+1)*10^(-9);
%     B_body_ctrl(:,i)=B_body_ctrl(:,i)+sqrt(R)*randn(size(B_body_ctrl(:,i)));
%     R_BO = R_OB';
% 
%     w_b_io = R_OB(:,2)*Const.w_o; 
%     w_b_ob(:,i) = x(5:7) - w_b_io; % Calculating angular rate of satellite relative to ECI frame
%     w_b_ob_magn(i) = norm(w_b_ob(:,i));
%     %% Bdot
%     if i==1
%         Bdot_body(:,i)=[0;0;0];
%         M=zeros(3,1);
%         Mag(:,i)= M;
%     elseif i>1
%         Bdot_body(:,i) = (B_body_ctrl(:,i) - B_body_ctrl(:,i-1))/dt;
% 
%         %% Torque Calculation
%         M = -Kp.*Bdot_body(:,i)./norm(B_body_ctrl(:,i));    % Calculate dipole moment to produce
%         M = mtq_scaling(M, mtq_max);                             % MTQ scaling 
%         Mag(:,i)= M;
%     end
%         T_magnetic = cross(M,B_body_ctrl(:,i));
%     for j=1:(dt/dt_model*0.3)
% 
%         q_orbit_eci=dcm2quat(Orbit2ECI_DCM(nodem(1,(i-1)*dt/dt_model+j), ...
%             inclm(1,(i-1)*dt/dt_model+j),argpm(1,(i-1)*dt/dt_model+j)+mm(1,(i-1)*dt/dt_model+j)));
%         q_eci_body = x(1:4);
%         q_orbit_body = quatProd(q_orbit_eci,q_eci_body);
%         R_OB = quat2dcm(q_orbit_body'); % Calculating the transformation matrix from orbit to body frame
%         R_BO = R_OB';
%         B_body=R_OB*mag_field_orbit(:,(i-1)*dt/dt_model+j)*10^(-9);
%         [T_disturbances, disturbancesMessage] = disturbances(R_BO, quat2eul(q_orbit_body'), ...
%             sun_pos_orbit(:,(i-1)*dt/dt_model+j), B_body, setDisturbances);
%         x = real_model.stateTransFun(x, stateTransCookieFinal(T_disturbances)); 
%     end
% 
%     for j=1:(dt/dt_model*0.7)
% 
%         q_orbit_eci=dcm2quat(Orbit2ECI_DCM(nodem(1,(i-1)*dt/dt_model+(dt/dt_model*0.3)+j), ...
%             inclm(1,(i-1)*dt/dt_model+(dt/dt_model*0.3)+j), ...
%             argpm(1,(i-1)*dt/dt_model+(dt/dt_model*0.3)+j)+mm(1,(i-1)*dt/dt_model+(dt/dt_model*0.3)+j)));
%         q_eci_body = x(1:4);
%         q_orbit_body = quatProd(q_orbit_eci,q_eci_body);
%         R_OB = quat2dcm(q_orbit_body'); % Calculating the transformation matrix from orbit to body frame
%         R_BO = R_OB';
%         B_body=R_OB*mag_field_orbit(:,(i-1)*dt/dt_model+(dt/dt_model*0.3)+j)*10^(-9);
%         [T_disturbances, disturbancesMessage] = disturbances(R_BO, quat2eul(q_orbit_body'), ...
%             sun_pos_orbit(:,(i-1)*dt/dt_model+(dt/dt_model*0.3)+j), ...
%             B_body, setDisturbances);
%         torq = T_magnetic + T_disturbances; 
%         x = real_model.stateTransFun(x, stateTransCookieFinal(torq)); 
%     end
%     t = t + dt;
% end

%% Simulation Parameters

fprintf("\n\n\nSIMULATION PARAMETERS\n");
fprintf("--------------------------------------------------------\n");
% fprintf("Initial Angular velocity: [%f, %f, %f]\n", w_init(1), w_init(2), w_init(3));
% fprintf("Initial Orientation: [%f, %f, %f]\n", angles_init(1), angles_init(2), angles_init(3));
fprintf("Detumbling Gain: %f\n", Kp); 
fprintf("External Torques active: %s\n", disturbancesMessage);
fprintf("Simulation ran for %f orbits\n", orbits);
fprintf("Maximum magnetic moment of MTQs set to: %f Am^2\n", Const.mtq_max);
fprintf("--------------------------------------------------------");
fprintf("\n\n\n");
fprintf("SIMULATION RESULTS\n");
fprintf("--------------------------------------------------------\n");
fprintf("Final Angular velocity: [%f, %f, %f]\n", w_b_ob(1,i), w_b_ob(2,i), w_b_ob(3,i));
fprintf("--------------------------------------------------------");
fprintf("\n\n\n");

%%  Plotting the Angular Velocities

    figure()
    subplot(3,1,1)
    plot(1:plotter_step:length(Time),w_b_ob(1,1:plotter_step:end))
    title('Angular Velocities','interpreter','latex', 'fontsize',17);
    %xlabel('Time [s]');
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    %ylabel('Velocity [rad/sec]');
    %ylabel(['\omega1 [rad/sec]'], 'interpreter','latex', 'fontsize',14);
    ylabel(['$\omega_' num2str(1) '$' '[rad/sec]'], 'interpreter','latex', 'fontsize',14);
    grid on
    subplot(3,1,2)
    plot(1:plotter_step:length(Time),w_b_ob(2,1:plotter_step:end))
    %title('Angular Velocities-y','interpreter','latex', 'fontsize',17);
    %xlabel('Time [s]');
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    %ylabel('Velocity [rad/sec]');
    %ylabel(['\omega2 [rad/sec]'], 'interpreter','latex', 'fontsize',14);
    ylabel(['$\omega_' num2str(2) '$' '[rad/sec]'], 'interpreter','latex', 'fontsize',14);
    grid on
    subplot(3,1,3)
    plot(1:plotter_step:length(Time),w_b_ob(3,1:plotter_step:end))
    %title('Angular Velocities-z','interpreter','latex', 'fontsize',17);
    %xlabel('Time [s]');
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    %ylabel(['\omega3 [rad/sec]'], 'interpreter','latex', 'fontsize',14);
    ylabel(['$\omega_' num2str(3) '$' '[rad/sec]'], 'interpreter','latex', 'fontsize',14);
    %ylabel('Velocity [rad/sec]');
    grid on

%% Plotting the angular velocity using Bdot 
    figure() 
    subplot(3,1,1) 
    plot(1:plotter_step:length(Time),w_b_ob_Bdot(1,1:plotter_step:end)) 
    title('Angular Velocity using Bdot') 
    xlabel('Time [s]'); 
    ylabel('Angular Velocity-x [rad/sec]'); 
    grid on; 
    subplot(3,1,2) 
    plot(1:plotter_step:length(Time),w_b_ob_Bdot(2,1:plotter_step:end)) 
    xlabel('Time [s]'); 
    ylabel('Angular Velocity-y [rad/sec]'); 
    grid on; 
    subplot(3,1,3) 
    plot(1:plotter_step:length(Time),w_b_ob_Bdot(3,1:plotter_step:end)) 
    xlabel('Time [s]'); 
    ylabel('Angular Velocity-z [rad/sec]'); 
    grid on; 

%% Bdot 

    figure()
    subplot(3,1,1)
    plot(1:plotter_step:length(Time),Bdot_body(1,1:plotter_step:end))
    title('Bdot-x');
    xlabel('Time[s]');
    ylabel('Bdot [T/sec)]');
    subplot(3,1,2)
    plot(1:plotter_step:length(Time),Bdot_body(2,1:plotter_step:end))
    title('Bdot-y');
    xlabel('Time[s]');
    ylabel('Bdot [T/sec)]');
    subplot(3,1,3)
    plot(1:plotter_step:length(Time),Bdot_body(3,1:plotter_step:end))
    title('Bdot-z');
    xlabel('Time[s]');
    ylabel('Bdot [T/sec)]');
%%  Plotting the Angular Velocity Magnitude

    figure()
    plot(1:plotter_step:length(Time),w_b_ob_magn(1:plotter_step:end))
    title('Angular Velocity Magnitude');
    xlabel('Time [s]');
    ylabel('|Velocity| [rad/sec]');
    
%%  Plotting the produced Torques

%     figure()
%     subplot(3,1,1)
%     plot(1:plotter_step:reps,tau(1,1:plotter_step:end))
%     title('Torques-x')
%     subplot(3,1,2)
%     plot(1:plotter_step:reps,tau(2,1:plotter_step:end))
%     title('Torques-y')
%     subplot(3,1,3)
%     plot(1:plotter_step:reps,tau(3,1:plotter_step:end))
%     title('Torques-z')
%     figure()
%     subplot(3,1,1)
%     plot(1:plotter_step:reps,tau(1,1:plotter_step:end)-T_disturbances(1,1:plotter_step:end))
%     title('Torques-x')
%     subplot(3,1,2)
%     plot(1:plotter_step:reps,tau(2,1:plotter_step:end)-T_disturbances(2,1:plotter_step:end))
%     title('Torques-y')
%     subplot(3,1,3)
%     plot(1:plotter_step:reps,tau(3,1:plotter_step:end)-T_disturbances(3,1:plotter_step:end))
%     title('Torques-z')
    
%%  Plotting the produced dipole moment

    figure()
    subplot(3,1,1)
    plot(1:plotter_step:length(Time),Mag(1,1:plotter_step:end))
    title('Magnetic Dipole-x');
    xlabel('Time [s]');
    ylabel({'Magnetic'; 'Dipole [Am^2]'});
    subplot(3,1,2)
    plot(1:plotter_step:length(Time),Mag(2,1:plotter_step:end))
    title('Magnetic Dipole-y');
    xlabel('Time [s]');
    ylabel({'Magnetic'; 'Dipole [Am^2]'});
    subplot(3,1,3)
    plot(1:plotter_step:length(Time),Mag(3,1:plotter_step:end))
    title('Magnetic Dipole-z');
    xlabel('Time [s]');
    ylabel({'Magnetic'; 'Dipole [Am^2]'});
    
%% Plotting Magnitudes Comparison [Disturbances, MTQs Torques]
%     figure()
%     plot(1:plotter_step:reps, T_disturbances_magn(1:plotter_step:end))
%     title('Magnitude of Disturbances')
%     ylabel('Torque [Nm]')
%     xlabel('Time [s]')
    
%     figure()
%     plot(1:plotter_step:reps, tau_magn(1:plotter_step:end))
%     title('Magnitude of MTQs Torques')
%     xlabel('Time [s]')
%     ylabel('Torque [Nm]')
    
%% Plotting disturbances

%     figure()
%     subplot(3,1,1)
%     plot(1:plotter_step:reps,T_disturbances(1,1:plotter_step:end))
%     ylabel('Torque [Nm]')
%     title('Disturbances-x')
%     subplot(3,1,2)
%     plot(1:plotter_step:reps,T_disturbances(2,1:plotter_step:end))
%     ylabel('Torque [Nm]')
%     title('Disturbances-y')
%     subplot(3,1,3)
%     plot(1:plotter_step:reps,T_disturbances(3,1:plotter_step:end))
%     title('Disturbances-z')
%     ylabel('Torque [Nm]')
