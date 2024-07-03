% ======================================================================== %
%   This function calculates the residual magnetic torque
% 
%   Inputs:
%     Area       - Satellite's area projected to the sun
%     cosines    - Cosine of angle between body frame axes and orbit frame z-axis
%     B_body     - Magnetic field expressed in body frame
% 
%   Ouputs:
%     tau_rm     - Residual magnetic torque
%     rm         - Residual magnetic moment
% 
% ======================================================================== %

function [tau_rm, rm] = residual_magnetic_torque(Area, cosines, B_body)

    % mean_rm_base = [0.05 0.05 0.05]'; 
    mean_rm_base = [0.01 0.01 0.01]';
    maximum_area = [0.034 0.034 0.01];
    
    sign_vector = sign(-cosines); 
    
    rm = mean_rm_base + (0.005*(Area./maximum_area).*sign_vector)' + 0.0005*rand(3,1);
    
    tau_rm = cross(rm, B_body);
end

