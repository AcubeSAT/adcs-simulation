% ======================================================================== %
%   This function calculates the disturbances that act upon the spacecraft.
% 
%   Inputs:
%     R_BO               - Transformation matrix from body to orbit frame
%     sun_vector         - Unit vector from satellite to sun expressed in orbit frame
%     B_body             - Magnetic field expressed in body frame
%     selector           - Message to select the active disturbances
% 
%   Ouputs:
%     T_dist             - Total disturbances torques
%     returnMessage      - Message that indicates which disturbance is active
% 
% ======================================================================== %

function [T_dist, returnMessage] = disturbances_bdot(R_BO, sun_vector, B_body, selector)

    Const = constants(); 
    R_OB = R_BO';
    
    %% Gravitational Torque

    [tau_g] = gravitational_torque(R_BO, Const.w_o, Const.I);

    %% Aerodynamic Drag

    [tau_ad] = aerodynamic_drag(R_BO, Const.p, Const.Cm, Const.v_satellite);
    
    %% Residual Magnetic Moment

    rm = [0.05 0.05 0.05]'; % An estimation of the residual moment of the satellite

    tau_rm = cross(rm, B_body);

    %% Solar Pressure

    [tau_sp, ~, ~] = solar_pressure(R_OB, sun_vector, Const.Cm);
    
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