% =========================================================================%
%   Function for calculating the voltage, amperage and power of the magnetorquers
%
%   Inputs:
%       m              - Magnetic Dipole provided to magnetorquers
%
%   Ouputs:
%       V_mtq          - Voltage on MTQs
%       I_mtq          - Amperage on MTQs
%       P_thermal_mtq  - Power consumed by MTQs
% =========================================================================%

function [V_mtq, I_mtq, P_thermal_mtq] = mtq_model(m)

    % Global constants may be taken as input from the function.

    global R_coils;
    global N_coils;
    global A_coils;
    global mt;
    % global K_inv

    % Array initialization.

    I_mtq = zeros(3, 1);
    V_mtq = zeros(3, 1);
    P_thermal_mtq = zeros(3, 1);

    %% V, I, W Calculation

    I_mtq(1, 1) = m(1) / (N_coils(1) * A_coils(1) * mt(1));
    I_mtq(2, 1) = m(2) / (N_coils(2) * A_coils(2) * mt(2));
    I_mtq(3, 1) = m(3) / (N_coils(3) * A_coils(3) * mt(3));

    V_mtq(1, 1) = I_mtq(1, 1) * R_coils(1);
    V_mtq(2, 1) = I_mtq(2, 1) * R_coils(2);
    V_mtq(3, 1) = I_mtq(3, 1) * R_coils(3);

    P_thermal_mtq(1, 1) = I_mtq(1, 1) * V_mtq(1, 1);
    P_thermal_mtq(2, 1) = I_mtq(2, 1) * V_mtq(2, 1);
    P_thermal_mtq(3, 1) = I_mtq(3, 1) * V_mtq(3, 1);

    % M = cross(B_body, T_desired)/norm(B_body,2)^2;
    % T_magnetic = cross(M, B_body);
     
    % V_mtq(:, i) = K_inv * cross(T_magnetic(:, i), B_body);
    % I_mtq(1, i) = V_mtq(1, i) / R_coils(1);
    % I_mtq(2, i) = V_mtq(2, i) / R_coils(2);
    % I_mtq(3, i) = V_mtq(3, i) / R_coils(3);

end