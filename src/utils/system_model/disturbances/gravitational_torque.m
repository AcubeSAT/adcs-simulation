% ======================================================================== %
%   This function calculates the gravitational torque
%
%   Inputs:
%     R_BO      - Rotation matrix from orbit to body frame
%     w_o       - Satellite angular velocity relative to Earth
%     I         - Inertia matrix
%
%   Ouputs:
%     tau_g     - Gravitational torque
%
% ======================================================================== %

function [tau_g] = gravitational_torque(R_BO, w_o, I)

    nadir = R_BO(:, 1);

    tau_g = 3 * w_o^2 * cross(nadir, I*nadir);

end