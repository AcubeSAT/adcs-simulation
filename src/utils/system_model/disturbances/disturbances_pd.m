% ======================================================================== %
%   This function calculates the disturbances that act upon the spacecraft.
%
%   Inputs:
%     q_orbit_body       - Quaternion that expresses the rotation from the orbit to the satellite body frame
%     sun_vector         - Unit vector from satellite to sun expressed in orbit frame
%     B_orbit            - Magnetic field expressed in orbit frame
%     selector           - Message to select the active disturbances
%
%   Ouputs:
%     T_dist             - Total disturbances torques
%     tau_g              - Gravitational torque
%     tau_ad             - Aerodynamic torque
%     tau_rm             - Residual magnetic torque
%     tau_sp             - Solar pressure torque
%     rm                 - Residual magnetic moment
%     returnMessage      - Message that indicates which disturbance is active
%
% ======================================================================== %
function [T_dist, rm, returnMessage, tau_ad, tau_rm, tau_sp, tau_g] = disturbances_pd(q_orbit_body, sun_vector, B_orbit, selector)

    Const = constants();

    R_OB = quat2dcm(q_orbit_body');
    R_BO = R_OB';
    B_body = R_OB * B_orbit;

    %% Gravitational Torque

    [tau_g] = gravitational_torque(R_BO, Const.w_o, Const.I);

    %% Aerodynamic Drag

    [tau_ad] = aerodynamic_drag(R_BO, Const.p, Const.Cm, Const.v_satellite);

    %% Solar Pressure

    [tau_sp, Area, cosines] = solar_pressure(R_OB, sun_vector, Const.Cm,Const.reflectance_factor);

    %% Residual Magnetic Moment

    [tau_rm, rm] = residual_magnetic_torque(Area, cosines, B_body);

    %% Total Torques

    T_total_dist = tau_g + tau_ad' + tau_rm + tau_sp;

    %% Return selector

    if (selector == "tau_g")

        T_dist = tau_g;
        returnMessage = "Only Gravity Gradient disturbance active.";

    elseif (selector == "tau_ad")

        T_dist = tau_ad;
        returnMessage = "Only Aerodynamic Torque active.";

    elseif (selector == "tau_rm")

        T_dist = tau_rm;
        returnMessage = "Only Residual Magnetic Moment disturbance active.";

    elseif (selector == "tau_sp")

        T_dist = tau_sp;
        returnMessage = "Only Solar Pressure disturbance active.";

    elseif (selector == "total")

        T_dist = T_total_dist;
        returnMessage = "All external torques active.";

    elseif (selector == "zero")

        T_dist = 0;
        returnMessage = "Disturbances set to zero.";

    end

end