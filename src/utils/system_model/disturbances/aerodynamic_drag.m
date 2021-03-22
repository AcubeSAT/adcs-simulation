% ======================================================================== %
%   This function calculates the aerodynamic drag that acts upon
%   the spacecraft.
% 
%   Inputs:
%     R_BO          - Transformation matrix from body to orbit frame
%     p             - Atmospheric density
%     Cm            - Center of mass
%     v_satellite   - Satellite's velocity in orbit
% 
%   Ouputs:
%     tau_ad        - Aerodynamic torque
% 
% ======================================================================== %

function [tau_ad] = aerodynamic_drag(R_BO, p, Cm, v_satellite)

    aerodynamic_constant = 2;

    z_ob = R_BO(:,3);                                   % Orbit frame z axis expressed in the body frame

    proj_Xb_Zo = ([1 0 0]*z_ob)*z_ob';                  % Projection of each body frame axis unit vector to orbit frame z axis unit vector
    proj_Yb_Zo = ([0 1 0]*z_ob)*z_ob';
    proj_Zb_Zo = ([0 0 1]*z_ob)*z_ob';

    proj_Xb_XYo = [1 0 0]-proj_Xb_Zo;                   % Projection of each body frame axis unit vector to orbit frame x,y plane
    proj_Yb_XYo = [0 1 0]-proj_Yb_Zo;
    proj_Zb_XYo = [0 0 1]-proj_Zb_Zo;

    Ax = 0.034*norm(cross(proj_Yb_XYo, proj_Zb_XYo));   % Surface projections to orbit frame x,y plane
    Ay = 0.034*norm(cross(proj_Xb_XYo, proj_Zb_XYo));
    Az = 0.01*norm(cross(proj_Xb_XYo, proj_Yb_XYo));

    uz_o = [0 0 1];

    cos_Xb_Xo = R_BO(1,:)*uz_o';                        % Cosine of angle between body frame axes and orbit frame z axis
    cos_Yb_Yo = R_BO(2,:)*uz_o';
    cos_Zb_Zo = R_BO(3,:)*uz_o';

    atmospheric_pressure_center = diag([0.05 0.05 0.17]);
    atmospheric_pressure_center(1,:) = sign(cos_Xb_Xo)*atmospheric_pressure_center(1,:); 
    atmospheric_pressure_center(2,:) = sign(cos_Yb_Yo)*atmospheric_pressure_center(2,:); 
    atmospheric_pressure_center(3,:) = sign(cos_Zb_Zo)*atmospheric_pressure_center(3,:); 

    T1 = 0.5 * p * aerodynamic_constant * Ax * v_satellite^2 * cross(z_ob', atmospheric_pressure_center(1,:)'-Cm);
    T2 = 0.5 * p * aerodynamic_constant * Ay * v_satellite^2 * cross(z_ob', atmospheric_pressure_center(2,:)'-Cm);
    T3 = 0.5 * p * aerodynamic_constant * Az * v_satellite^2 * cross(z_ob', atmospheric_pressure_center(3,:)'-Cm);

    tau_ad = T1+T2+T3;

end

