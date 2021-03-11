% *************************************************************************
% @file           : disturbances.m
% @brief          : Function that calculates the external torques that act
%                   on the satellite.
% *************************************************************************
% @attention
% 
% <h2><center>&copy; Copyright (c) Aristotle Space & Aeronautics Team.
% All rights reserved.</center></h2>
%  
% This project is licensed under the GNU General Public License v3.0.
%  
% *************************************************************************

function [T_dist, rm, returnMessage] = disturbances_pd(q_orbit_body, sun_vector, B, selector)

Const = constants(); % Loading some constant parameters.

R_OB = quat2dcm(q_orbit_body');
angles = quat2eul(q_orbit_body');
R_BO = R_OB';
B = R_OB*B; % orbit to body

%% Gravitational Torque

z_0 = R_BO(:,1); % Unit vector towards nadir

tau_g = 3 * Const.w_o^2 * cross(z_0, Const.I * z_0);

%% Aerodynamic Drag

Cd = 2;          % Aerodynamic Constant
p = Const.p;     % Atmospheric density
Cm = Const.Cm;    % Center of mass
% A_s = 0.01*cos(angles(2))*cos(angles(3)) + 0.034*sin(angles(2))*cos(angles(3)) + 0.034*sin(angles(3)); % Drag area, calculated by Trajectory

x_0 = R_BO(:,3);  % Velocity vector

projxb_zo = ([1 0 0]*x_0)*x_0'; % Projection of each body frame axis unit vector to orbit frame z-axis unit vector
projyb_zo = ([0 1 0]*x_0)*x_0';
projzb_zo = ([0 0 1]*x_0)*x_0';

projXb_XYo = [1 0 0]-projxb_zo; % Projection of each body frame axis unit vector to orbit frame x,y plane
projYb_XYo = [0 1 0]-projyb_zo;
projZb_XYo = [0 0 1]-projzb_zo;

Ax = 0.034*norm(cross(projYb_XYo,projZb_XYo)); % Surface projections to orbit frame x,y plane
Ay = 0.034*norm(cross(projXb_XYo,projZb_XYo));
Az = 0.01*norm(cross(projXb_XYo,projYb_XYo));

uz_o = [0 0 1];

cos_Xb_Xo = R_BO(1,:)*uz_o'; % Cosine of angle between Body frame axes and Orbit frame z-axis
cos_Yb_Yo = R_BO(2,:)*uz_o';
cos_Zb_Zo = R_BO(3,:)*uz_o';

Cpa = diag([0.05 0.05 0.17]);     % Center of atmospheric pressure for each side

Cpa(1,:) = sign(cos_Xb_Xo)*Cpa(1,:); 
Cpa(2,:) = sign(cos_Yb_Yo)*Cpa(2,:); 
Cpa(3,:) = sign(cos_Zb_Zo)*Cpa(3,:); 

T1 = 0.5 * p * Cd * Ax * Const.v_satellite^2 * cross(x_0', Cpa(1,:)'-Cm);
T2 = 0.5 * p * Cd * Ay * Const.v_satellite^2 * cross(x_0', Cpa(2,:)'-Cm);
T3 = 0.5 * p * Cd * Az * Const.v_satellite^2 * cross(x_0', Cpa(3,:)'-Cm);

tau_ad = T1+T2+T3;



%% Solar Pressure

Fs = 1367;  % Solar Constant [W/m^2]
c  = 3e8;   % Speed of light [m/s]
q  = 0.6;   % Reflectance factor
y_0 = -R_OB*sun_vector; % Sun vector in body frame

y_0 = y_0 /norm(y_0);

projxb_zo = ([1 0 0]*y_0)*y_0';
projyb_zo = ([0 1 0]*y_0)*y_0';
projzb_zo = ([0 0 1]*y_0)*y_0';

projXb_XYo = [1 0 0]-projxb_zo; % Projection of each body frame axis unit vector to orbit frame x,y plane
projYb_XYo = [0 1 0]-projyb_zo;
projZb_XYo = [0 0 1]-projzb_zo;

Ax = 0.034*norm(cross(projYb_XYo,projZb_XYo)); % Surface projections to orbit frame x,y plane
Ay = 0.034*norm(cross(projXb_XYo,projZb_XYo));
Az = 0.01*norm(cross(projXb_XYo,projYb_XYo));


cos_Xb_Xo = [1 0 0]*y_0;
cos_Yb_Yo = [0 1 0]*y_0;
cos_Zb_Zo = [0 0 1]*y_0;

Csp = diag([0.05 0.05 0.17]); % Center of solar pressure 

Csp(1,:) = sign(-cos_Xb_Xo)*Csp(1,:); 
Csp(2,:) = sign(-cos_Yb_Yo)*Csp(2,:); 
Csp(3,:) = sign(-cos_Zb_Zo)*Csp(3,:); 

T1 = (Fs/c)*Ax*(1+q)*cross(y_0,Csp(1,:)'-Cm);
T2 = (Fs/c)*Ay*(1+q)*cross(y_0,Csp(2,:)'-Cm);
T3 = (Fs/c)*Az*(1+q)*cross(y_0,Csp(3,:)'-Cm);

tau_sp = (T1+T2+T3);


%% Residual Magnetic Moment

rm_base = [0.05 0.05 0.05]'; % Mean value
A = [Ax Ay Az];
sign_vector = [sign(-cos_Xb_Xo) sign(-cos_Yb_Yo) sign(-cos_Zb_Zo)]; 
maximum_area = [0.034 0.034 0.01];
rm = rm_base + (0.005*(A./maximum_area).*sign_vector)' + 0.0005*rand(3,1);
tau_rm = cross(rm, B);

%% Total Torques

T_total_dist = tau_g + tau_ad' + tau_rm + tau_sp;

%% Return selector

if(selector == "tau_g")
    
    T_dist = tau_g;
    returnMessage = "Only Gravity Gradient disturbance active.";

elseif(selector == "tau_ad")
    
    T_dist = tau_ad;
    returnMessage = "Only Aerodynamic Torque active.";  
    
elseif(selector == "tau_rm")
    
    T_dist = tau_rm;
    returnMessage = "Only Residual Magnetic Moment disturbance active.";
    
elseif(selector == "tau_sp")
    
    T_dist = tau_sp;
    returnMessage = "Only Solar Pressure disturbance active.";
    
elseif(selector == "total")
    
    T_dist = T_total_dist;
    returnMessage = "All external torques active.";
    
elseif(selector == "zero")    

    T_dist = 0;
    returnMessage = "Disturbances set to zero.";
    
end

end