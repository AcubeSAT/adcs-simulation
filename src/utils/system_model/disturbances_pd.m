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

function [T_dist, returnMessage] = disturbances_pd(q_orbit_body, sun_vector, B, selector)

Const = constants(); % Loading some constant parameters.

R_OB = quat2dcm(q_orbit_body');
angles = quat2eul(q_orbit_body');
R_BO = R_OB';
B = R_OB*B; % orbit to body

%% Gravitational Torque

z_0 = R_BO(:, 1); % Unit vector towards nadir

tau_g = 3 * Const.w_o^2 * cross(z_0, Const.I * z_0);

%% Aerodynamic Drag

Cd  = 2;          % Aerodynamic Constant
p   = 10e-13;     % Atmospheric density
A_s = 0.01*cos(angles(2))*cos(angles(3)) + 0.034*sin(angles(2))*cos(angles(3)) + 0.034*sin(angles(3)); % Drag area, calculated by Trajectory
x_0 = R_BO(:,3);  % Velocity vector
Cpa = [0 0 0.17]';

tau_ad = 0.5 * p * Cd * A_s * Const.v_satellite^2 * cross(x_0, Cpa);

%% Residual Magnetic Moment

rm = [0.05 0.05 0.05]'; % An estimation of the residual moment of the satellite

tau_rm = cross(rm, B);

%% Solar Pressure

Fs = 1367;  % Solar Constant [W/m^2]
c  = 3e8;   % Speed of light [m/s]
Ap = 0.034; % Projected Area [m^2]
q  = 0.6;   % Reflectance factor
angle_inc_cos = angle_calc(angles, sun_vector);
x_0 = R_BO(:,3);  % Velocity vector
Csp = [0 0.05 0]'; % Center of solar pressure 

tau_sp = (Fs/c)*Ap*(1+q)*angle_inc_cos*cross(x_0,Csp);

%% Total Torques

T_total_dist = tau_g + tau_ad + tau_rm + tau_sp;

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