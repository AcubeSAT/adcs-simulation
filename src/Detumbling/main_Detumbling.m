% ========================================================================
%   Main function for Detumbling Mode simulation.
% ========================================================================

close all;
clc;
clear variables;

%% Initialize Parameters Script
Const = constants();
setParamsFinal_Detumbling();
Time = 0:cycle_duration:Total_simulation_time;
x_data = zeros(7, length(Time)); % Real state
x = x0;
t = 0;
B_body = zeros(3, length(Time));                    % Earth magnetic field satellite in body frame
w_b_ob = zeros(3, length(Time));                    % Angular velocity from orbit to body frame expressed in body frame
Bdot_body = zeros(3, length(Time));                 % Bdot estimation in body frane
w_b_ob_magn = zeros(1, length(Time));               % Angular velocity magnitude
Mag = zeros(3, length(Time));                       % Commanded magnetic dipole moment
w_b_ob_Bdot = zeros(3, length(Time));               % Angular velocity estimation using Bdot
q_orbit_body_data = zeros(4, length(Time)); 
threshold_times = 0;
threshold_exceptions = 0;
nominal_activation_matrix = zeros(2, length(Time)); % 1 - Bdot, 2 - nonBdot


global R_coils; R_coils = Const.R_coils;
global N_coils; N_coils = Const.N_coils;
global A_coils; A_coils = Const.A_coils;
global mt; mt = Const.mt;

for current_cycle = 1:length(Time) %Main loop

    %% Angular velocity and Magnetic field expressed in body frame
    q_orbit_eci = dcm2quat(Orbit2ECI_DCM(nodem(1, current_cycle), inclm(1, current_cycle), argpm(1, current_cycle)+mm(1, current_cycle)));
    q_eci_body = x(1:4);
    q_orbit_body = quatProd(q_orbit_eci, q_eci_body);
    q_orbit_body_data(:, current_cycle) = q_orbit_body;
    R_OB = quat2dcm(q_orbit_body'); % Calculating the transformation matrix from orbit to body frame
    B_body(:, current_cycle) = R_OB * mag_field_orbit(:, current_cycle) * 10^(-9);
    B_body(:, current_cycle) = B_body(:, current_cycle) + sqrt(R) * randn(size(B_body(:, current_cycle))); % Adding white noise to magnetometer measurements
    R_BO = R_OB';
    w_b_io = R_OB(:, 3) * Const.w_o;
    w_b_ob(:, current_cycle) = x(5:7) - w_b_io; % Calculating angular rate of satellite relative to ECI frame
    w_b_ob_magn(current_cycle) = norm(w_b_ob(:, current_cycle));

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
    B_body_2 = B_body_2 + sqrt(R) * randn(size(B_body(:, current_cycle))); % Second measurment from magnetometer

    %% Bdot calculation
    Bdot_body(:, current_cycle) = (B_body_2 - B_body(:, current_cycle)) / 0.1;

    %% Torque Calculation
    M = -Kp .* Bdot_body(:, current_cycle) ./ norm(B_body(:, current_cycle)); % Calculate dipole moment to produce
    M = mtq_scaling(M, Const.mtq_max); % MTQ scaling
    Mag(:, current_cycle) = M;
    T_magnetic = cross(M, B_body(:, current_cycle));

    %% Calculate V, I, P of MTQ's
    [V_mtq, I_mtq, P_thermal_mtq] = mtq_model(M);

    %% Estimation of angular velocity using Bdot
    w_b_ob_Bdot(:, current_cycle) = skew(B_body(:, current_cycle)) * (-Bdot_body(:, current_cycle) ...
        -B_body(:, current_cycle)) / (B_body(:, current_cycle)' * B_body(:, current_cycle));

    %% Calculate when Nominal is ready to be activated
    if current_cycle > 1
        [bdot_activation, nonBdot_activation, threshold_times, threshold_exceptions] = trigger_D2N(threshold_times, threshold_exceptions, w_b_ob(:, current_cycle), w_b_ob_Bdot(:, current_cycle), w_b_ob_Bdot(:, current_cycle-1),Const.D2N_threshold,total_limit,exceptions_limit);
        nominal_activation_matrix(1, current_cycle) = bdot_activation;
        nominal_activation_matrix(2, current_cycle) = nonBdot_activation;
    end

    %%  Disturbances for 2 timesteps to acquire magnetometer measurements
    [T_disturbances, ~] = disturbances_bdot(R_BO, sun_pos_orbit(:, current_cycle), B_body(:, current_cycle), setDisturbances);
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
title('Angular Velocities', 'interpreter', 'latex', 'fontsize', 17);
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

%% Plotting the angular velocity using Bdot
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

%% Bdot
figure()
subplot(3, 1, 1)
plot(1:plotter_step:length(Time), Bdot_body(1, 1:plotter_step:end))
title('Bdot-x');
xlabel('Time[s]');
ylabel('Bdot [T/sec)]');
subplot(3, 1, 2)
plot(1:plotter_step:length(Time), Bdot_body(2, 1:plotter_step:end))
title('Bdot-y');
xlabel('Time[s]');
ylabel('Bdot [T/sec)]');
subplot(3, 1, 3)
plot(1:plotter_step:length(Time), Bdot_body(3, 1:plotter_step:end))
title('Bdot-z');
xlabel('Time[s]');
ylabel('Bdot [T/sec)]');

%%  Plotting the Angular Velocity Magnitude
figure()
plot(1:plotter_step:length(Time), w_b_ob_magn(1:plotter_step:end))
title('Angular Velocity Magnitude');
xlabel('Time [s]');
ylabel('|Velocity| [rad/sec]');

%%  Plotting the produced dipole moment
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

%% Plotting Nominal Activation Matrix
figure();
for i = 1:2
    subplot(2, 1, i);
    hold on;
    plot(Time(1:length(Time)), nominal_activation_matrix(i, 1:length(Time)), 'LineWidth', 1.5, 'Color', 'blue');
    if (i == 1), title('Nominal Activation', 'interpreter', 'latex', 'fontsize', 17); end
    if (i == 1), ylabel('Data after process', 'interpreter', 'latex', 'fontsize', 14); end
    if (i == 2), ylabel('Raw Data', 'interpreter', 'latex', 'fontsize', 14); end
    xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
    hold off;
    grid on;
end

%% Plotting eclipse
figure()
plot(1:length(eclipse), eclipse, 'LineWidth', 2.0, 'Color', 'blue');
title('Eclipse', 'interpreter', 'latex', 'fontsize', 17)
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
ylabel('Eclipse', 'interpreter', 'latex', 'fontsize', 14);
grid on;
if (i == 1), title('Umbral, Penumbral or no Eclipse', 'interpreter', 'latex', 'fontsize', 17); end

%% 3D plot
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
% 
% 
% face_colors = [0 1 1;  % Cyan for -Z
%                1 0 1;  % Magenta for +Z
%                0 1 1;  % Cyan for bottom (-Y)
%                1 0 1;  % Magenta for +X
%                0 1 1;  % Cyan for top (+Y)
%                0 1 1]; % Cyan for -X
% 
% 
% 
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
%  rect_prism = patch('Vertices', vertices, 'Faces', faces, ...
%                    'FaceVertexCData', face_colors, 'FaceColor', 'flat', ...
%                    'EdgeColor', 'black');
% 
% 
% % Initialize text handle for the current time and frame
% time_text = text(0, 5, 5, '', 'FontSize', 12, 'Color', 'black');  % Position the text above the plot area
% 
% % Initialize quiver objects for the vectors
%  B_quiver = quiver3(0, 0, 0, 0, 0, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Magnetic field vector (Blue)
% 
%  % Add legend
%  legend({'Satellite', 'Magnetic Field'}, 'Location', 'best');
% 
% % Animation loop
% for t = 1:length(Time)
% 
%     % Rotation matrix
%     q_orbit_body = q_orbit_body_data(:,t);
%     R_OB = quat2dcm(q_orbit_body');
% 
% 
%     % Apply rotation to each vertex
%     rotated_vertices = (R_OB' * vertices')';  % Transform vertices using R
% 
%     % Update the rectangular prism's vertices
%     set(rect_prism, 'Vertices', rotated_vertices);
% 
% 
%      %Update Magnetic Field Vector
% 
%     B_orbit =R_OB'*B_body(:,t)*1e5;   % Transform to orbit frame
%     
% 
%     set(B_quiver, 'XData', 0, 'YData', 0, 'ZData', 0, ...
%                   'UData', B_orbit(1), 'VData', B_orbit(2), 'WData', B_orbit(3));
% 
% 
%     % Update the time and frame label
%     set(time_text, 'String', sprintf('Time: %d', t));
% 
%     % Update the plot
%     drawnow;
% 
%     pause(0.05);
% end 
