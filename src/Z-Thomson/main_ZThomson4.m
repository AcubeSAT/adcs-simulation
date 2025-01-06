% ========================================================================
%   Main function for Z Thomson simulation.
% ========================================================================

close all;
clear;
clc;


%% Initialize Parameters Script

Const=constants();
Param = setParams_ZThomson(Const.I);
dt = Param.dt;
orbits = Param.orbits;
tf = Param.tf;
x0 = Param.x0;
x0_hat = Param.x0_hat;
real_model = Param.real_model;
%model = Param.model;

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
N_Timesteps= Param.N_Timesteps;


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




%% Initialize simulation parameters

Time = 0:dt:tf;
x_real = zeros(7,length(Time));             % Real state
q_ob_data = zeros(4,length(Time));
x = x0(1:7);
x_real(:,1) = x0(1:7);
t = 0;
number_of_cycles = floor(length(Time)/N_Timesteps);  % Number of cycles
timeflag_dz = 0;
init_AngVel_dz = 0;
init_accel_dz = 0;
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
bias_data = zeros(3,length(Time));
gyro_noise_data = zeros(3,length(Time));
Bbody_data = zeros(3,length(Time));
bdot_activation_matrix = zeros(2, length(Time));
theta_deg_arr_x = zeros(1,length(Time));
theta_deg_arr_y = zeros(1,length(Time));
theta_deg_arr_z = zeros(1,length(Time));
y_noise_data = zeros(9,length(Time));

real_bias=0;

%% Main continuous loop

for cycle_index = 1:number_of_cycles

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

        q_orbit_eci = dcm2quat(Orbit2ECI_DCM(Nodem, Inclm, Argpm+Mm));
        q_eci_body = x(1:4);
        q_orbit_body = quatProd(q_orbit_eci, q_eci_body);
        R_OB = quat2dcm(q_orbit_body');

     

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
        theta_deg_arr_x(current_timestep) = theta_deg_x;
        theta_deg_arr_y(current_timestep) = theta_deg_y;
        theta_deg_arr_z(current_timestep) = theta_deg_z;


        %% Propagate the system
        q_ob = quat_EB2OB(x(1:4), Nodem,Inclm,Argpm,Mm );
        [T_dist, ~,~,ad,r,sp,g] = disturbances_pd(q_ob,Sun_pos_orbit, Mag_field_orbit, setDisturbances);

        torq = T_dist;
        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,0,[0;0;0]));

        %% Matrices update

        tau_ad(:,current_timestep) = ad;
        tau_rm(:,current_timestep) = r;
        tau_sp(:,current_timestep) = sp;
        tau_g(:,current_timestep) = g;
        tau_dist(:,current_timestep) = T_dist;
        x_real(:,current_timestep)=x;
        Bbody_data(:,current_timestep) = y_real(1:3)*norm(mag_field_orbit(:,current_timestep)*10^(-9));
        bias_data(:,current_timestep) = real_bias;
        gyro_noise_data(:,current_timestep) = gyro_noise;
        q_ob_data(:,current_timestep) = q_ob;
        y_noise_data(:,current_timestep) = y_noise;

       

    end
       
    w_b_io = R_OB(:, 3) * Const.w_o;
    w_b_ob = x(5:7) - w_b_io; % Calculating angular rate of satellite relative to ECI frame
    w_b_ib =  w_b_ob + w_b_io;
   
    B_body_2 = y_noise_data(1:3,current_timestep-1);
    B_body = y_noise_data(1:3,current_timestep-2);

    Bz=acos(B_body(3)/norm(B_body));
    Bz2=acos(B_body_2(3)/norm(B_body_2));

    Bdot_body_z = (Bz2 - Bz) / 0.1;


    %% Thomson spin
    M(3,1)= Param.Kd * Bdot_body_z;
    if (abs(B_body(2))>abs(B_body(1)))
        M(1,1)= -Param.Ks*(abs(w_b_ib(3))-Param.w_ref)*sign(B_body(2));
        M(2,1)= 0;
    elseif (abs(B_body(2))<abs(B_body(1)))
        M(1,1)= 0;
        M(2,1)=Param.Ks*(abs(w_b_ib(3))-Param.w_ref)*sign(B_body(1));
    end  


     M = mtq_scaling(M, Const.mtq_max); % MTQ scaling
     T_magnetic = cross(M, B_body);


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

        q_orbit_eci = dcm2quat(Orbit2ECI_DCM(nodem(1, current_timestep), inclm(1, current_timestep), argpm(1, current_timestep)+mm(1, current_timestep)));
        q_eci_body = x(1:4);
        q_orbit_body = quatProd(q_orbit_eci, q_eci_body);
        R_OB = quat2dcm(q_orbit_body');

        %% Sensor Measurements

        y_real = real_model.msrFun(x,msrCookieFinal(Mag_field_eci,Sun_pos_eci,Eclipse,[0;0;0]));

        [gyro_noise,real_bias] = gyro_noise_func(real_bias,dt,sigma_u,sigma_v);
        y_noise(4:6) = y_real(4:6) + gyro_noise;

        
 
      
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
        theta_deg_arr_x(current_timestep) = theta_deg_x;
        theta_deg_arr_y(current_timestep) = theta_deg_y;
        theta_deg_arr_z(current_timestep) = theta_deg_z;

        %% Propagate the system

        q_ob = quat_EB2OB(x(1:4),Nodem,Inclm,Argpm,Mm);
        [T_dist, ~,~,ad,r,sp,g] = disturbances_pd(q_ob,Sun_pos_orbit, Mag_field_orbit, setDisturbances);

        torq = T_magnetic + T_dist;

        x = real_model.stateTransFun(x, stateTransCookieFinalNominal(torq,0,[0;0;0]));

        %% Matrices update

        tau_mtq(:,current_timestep) = T_magnetic;
        tau_ad(:,current_timestep) = ad;
        tau_rm(:,current_timestep) = r;
        tau_sp(:,current_timestep) = sp;
        tau_g(:,current_timestep) = g;
        tau_dist(:,current_timestep) = T_dist;
        x_real(:,current_timestep)=x;
        M_data(:,current_timestep) = M;
        bias_data(:,current_timestep) = real_bias;
        gyro_noise_data(:,current_timestep) = gyro_noise;
        Bbody_data(:,current_timestep) = y_real(1:3)*norm(mag_field_orbit(:,current_timestep)*10^(-9));
        q_ob_data(:,current_timestep) = q_ob;

        

    end
end




%% =============================== Errors and plots =============================================== %%

%% Eclipse plot

figure()
plot(1:length(eclipse),eclipse, 'LineWidth',2.0, 'Color','blue');
title('Eclipse', 'interpreter','latex', 'fontsize',17)
xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
ylabel('Eclipse', 'interpreter','latex', 'fontsize',14);
grid on;
if (i==1), title('Umbral, Penumbral or no Eclipse', 'interpreter','latex', 'fontsize',17);end



%%
figure();
for i=1:3
    subplot(3,1,i);
    plot(Time(1:end-1),x_real(4+i,1:length(Time)-1),'LineWidth',2.0, 'Color','blue');
    xlabel('Time [s]', 'interpreter','latex', 'fontsize',12);
    if (i==1), ylabel('$\omega_1$ [rad/sec]', 'interpreter','latex', 'fontsize',14); end
    if (i==2), ylabel('$\omega_2$ [rad/sec]', 'interpreter','latex', 'fontsize',14); end
    if (i==3), ylabel('$\omega_3$ [rad/sec]', 'interpreter','latex', 'fontsize',14); end
    ylabel(['$\omega_' num2str(i) '$' '[rad/sec]'], 'interpreter','latex', 'fontsize',14);
    if (i==1), legend('Angular Velocity');end
    if (i==1), title('Angular Velocities', 'interpreter','latex', 'fontsize',17); end
    grid on;
end


%%  Plotting the produced Torques

total_torques = tau_mtq;

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

%%
mtq_max_vector = Const.mtq_max';
mtq_max_vector = mtq_max_vector * ones(1,length(Bbody_data));
max_mtq_torque = abs(cross(mtq_max_vector,Bbody_data));
figure()
for i=1:3
    subplot(3,3,i);
    hold on;
    plot(Time(1:length(Time)), max_mtq_torque(i, 1:length(Time)), 'LineWidth',1.5, 'Color','blue');
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
    plot(Time(1:length(Time)), total_torques(i,1:length(Time)), 'LineWidth',1.5, 'Color','blue');
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
    plot(Time(22:length(Time)),tau_dist(i,22:length(Time)), 'LineWidth',1.5, 'Color','blue')
    if (i==2), title('Total disturbance torque', 'interpreter','latex', 'fontsize',17);end
    if (i==1), ylabel('X-axis [deg]'); end
    if (i==2), ylabel('Y-axis [deg]'); end
    if (i==3), ylabel('Z-axis [deg]'); end
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    hold off;
    grid on;
end

%% Max actuator Torque - Actuator Torque - Total Disturbances


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


%%  Plotting the Disturbances

figure()
subplot(3,1,1)
plot(Time(1:end-1),tau_dist(1,1:length(Time)-1), 'LineWidth',1.5, 'Color','blue')
title('Disturbances', 'interpreter','latex', 'fontsize',17)
ylabel('Disturbances-x [Nm]', 'interpreter','latex', 'fontsize',14)
xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
grid on;
subplot(3,1,2)
plot(Time(1:end-1),tau_dist(2,1:length(Time)-1), 'LineWidth',1.5, 'Color','blue')
ylabel('Disturbances-y [Nm]', 'interpreter','latex', 'fontsize',14)
xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
grid on;
subplot(3,1,3)
plot(Time(1:end-1),tau_dist(3,1:length(Time)-1), 'LineWidth',1.5, 'Color','blue')
ylabel('Disturbances-z [Nm]', 'interpreter','latex', 'fontsize',14)
xlabel('Time [s]', 'interpreter','latex', 'fontsize',12)
grid on;


    


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
