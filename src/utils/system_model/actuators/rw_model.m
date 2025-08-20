% =========================================================================%
%   Function for calculating the current state of the Reaction Wheel
%
%   Inputs:
%       T_rw_desired        - Desired torque
%       AngVel_rw           - Current angular velocity of the RW
%
%   Outputs:
%       V_rw                - Voltage of RW
%       I_rw                - Amperage of RW
%       P_thermal_rw        - Total power consumed on RW
%       AngVel_rw_new       - Angular velocity of RW on next timestep
%       acceleration        - Acceleration of RW
%       T_rw_total          - Total RW torque (with frictions)
%
% =========================================================================%

function [V_rw, I_rw, P_thermal_rw, AngVel_rw_new, acceleration, T_rw_total] = rw_model(T_rw_desired, AngVel_rw)
    % The angular velocity of the Reaction Wheel is given by a sensor
    % placed on the RW.

    timestep = 0.1;

    global Jw;
    global Km;
    global Kv;
    global b_friction;
    global c_friction;
    global Rb;

    % global Ai;
    % global T_friction;
    % global Kv;

    %% V, I, W Calculation
    
    T_rw_total = T_rw_desired; 
    
    acceleration = (T_rw_total)/Jw;      %In rad/sec^2
    
    I_rw = T_rw_desired / Km + AngVel_rw / (Kv * Rb);

    % % If T_friction is considered = 0, may be omitted
    % V_rw = Kv * AngVel_rw - Rb * (1/Km) * (Ai * (-T_desired(1) - T_friction));

    % % Initial approach of functionality in deadzone
    % f = I_rw * Km / c;
    %
    % if abs(f) < 1
    %     sat = f;
    % elseif f > 1
    %     sat = 1;
    % else
    %     sat = -1;
    % end
    %
    % if abs(AngVel_rw) < lim_dz
    %     if abs(sat) == 1
    %         I_rw = (Jw * acceleration + c * sat)/Km;
    %     else
    %         I_rw = 0;
    %         acceleration = 0;
    %     end
    % end

    AngVel_rw_new = AngVel_rw + acceleration * timestep; %In rad/s = rad/s + (rad/s^2)*sec
    acceleration = acceleration * 30 / pi; %In rpm/sec
    V_rw = I_rw * Rb;
    P_thermal_rw = I_rw * V_rw;

    % If the values need to be ploted outside the main loop of the controler
    %
    % if plots_on == 1
    %     %%  Plotting the Voltage of RW
    %
    %     figure()
    %     plot(t(1:10:end),V_rw(1:10:end))
    %     title('Voltage of RW')
    %
    %     %%  Plotting the Current of RW
    %
    %     figure()
    %     plot(t(1:10:end),I_rw(1:10:end))
    %     title('Current of RW')
    %
    %     %%  Plotting the Power of RW
    %
    %     figure()
    %     plot(t(1:10:end),W_rw(1:10:end))
    %     title('Power of RW')
    % end
end
