%% ========================================================================%%
%   Function for desaturating the RW in case that the angular velocity of
%   the RW exceeds a specified value.
%   In case this happens, part of the RW Torque is added on the Magnetic
%   Torque (if possible).
%
%   Inputs: Current Magnetic Torque 
%           RW Torque
%           RW Acceleration
%           RW Angular Velocity
%           Magnetic Field expressed on Body frame
%           Maximum Dipole provided by each Magnetorquer
%   Ouputs: Updated Magnetic Torque
%           Updated RW Torque
%
% =========================================================================%

function [T_mtq_effective, T_rw] = rw_saturation(T_magnetic_effective, T_rw, accel_rw, AngVel_rw, B_body, mtq_max)

    % The angular velocity of the Reaction Wheel is given by a sensor
    % placed on the RW.
    
    global A;
    global Jw;
    
    AngVel_rw_lim = 10000;
    T_mtq_effective = T_magnetic_effective;

    if abs(AngVel_rw) > AngVel_rw_lim && abs(T_rw(3))> 0
        T_added = [0; 0; A * Jw * accel_rw];
        T_magnetic = T_magnetic_effective + T_added;
        M = -cross(T_magnetic,B_body)/(norm(B_body))^2;
        M = mtq_scaling(M, mtq_max);
        T_mtq_effective = cross(M,B_body);
        
        T_rw = T_rw - (T_mtq_effective - T_magnetic_effective);
    end
end

