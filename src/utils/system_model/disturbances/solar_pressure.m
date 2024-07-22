% ======================================================================== %
%   This function calculates the solar pressure that acts upon
%   the spacecraft.
%
%   Inputs:
%     R_OB          - Transformation matrix from orbit to body frame
%     sun_vector    - Unit vector from satellite to sun expressed in orbit frame
%     Cm            - Center of mass
%     reflectance_factor - Reflectance factor of the spacecraft surface materials
%
%   Ouputs:
%     tau_ad        - Solar pressure torque
%     Area          - Satellite's area projected to the sun
%     cosines       - Cosine of angle between Body frame axes and orbit frame z-axis
%
% ======================================================================== %

function [tau_sp, Area, cosines] = solar_pressure(R_OB, sun_vector, Cm,reflectance_factor)

    Fs = 1367; % Solar Constant [W/m^2]
    c = 3e8; % Speed of light [m/s]
    

    y_0 = -R_OB * sun_vector; % Sun vector in body frame
    y_0 = y_0 / norm(y_0);

    proj_Xb_Zo = ([1, 0, 0] * y_0) * y_0';
    proj_Yb_Zo = ([0, 1, 0] * y_0) * y_0';
    proj_Zb_Zo = ([0, 0, 1] * y_0) * y_0';

    proj_Xb_XYo = [1, 0, 0] - proj_Xb_Zo; % Projection of each body frame axis unit vector to orbit frame x,y plane
    proj_Yb_XYo = [0, 1, 0] - proj_Yb_Zo;
    proj_Zb_XYo = [0, 0, 1] - proj_Zb_Zo;

    Ax = 0.034 * norm(cross(proj_Yb_XYo, proj_Zb_XYo)); % Surface projections to orbit frame x,y plane
    Ay = 0.034 * norm(cross(proj_Xb_XYo, proj_Zb_XYo));
    Az = 0.01 * norm(cross(proj_Xb_XYo, proj_Yb_XYo));
    Area = [Ax, Ay, Az];

    cos_Xb_Xo = [1, 0, 0] * y_0;
    cos_Yb_Yo = [0, 1, 0] * y_0;
    cos_Zb_Zo = [0, 0, 1] * y_0;
    cosines = [cos_Xb_Xo, cos_Yb_Yo, cos_Zb_Zo];

    solar_pressure_center = diag([0.05, 0.05, 0.17]);

    solar_pressure_center(1, :) = sign(-cos_Xb_Xo) * solar_pressure_center(1, :);
    solar_pressure_center(2, :) = sign(-cos_Yb_Yo) * solar_pressure_center(2, :);
    solar_pressure_center(3, :) = sign(-cos_Zb_Zo) * solar_pressure_center(3, :);

    T1 = (Fs / c) * Ax * (1 + reflectance_factor) * cross(y_0, solar_pressure_center(1, :)'-Cm);
    T2 = (Fs / c) * Ay * (1 + reflectance_factor) * cross(y_0, solar_pressure_center(2, :)'-Cm);
    T3 = (Fs / c) * Az * (1 + reflectance_factor) * cross(y_0, solar_pressure_center(3, :)'-Cm);

    tau_sp = (T1 + T2 + T3);
    
end
