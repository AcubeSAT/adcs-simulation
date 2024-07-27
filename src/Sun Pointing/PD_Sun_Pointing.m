% ======================================================================== %
%   Implementation of the PD controller that calculates the desired torque
%   to be applied in order to achieve sun pointing.
%
%   Inputs:
%     q_desired              - Desired quaternion
%     q_eci_body             - Quaternion that expresses the rotation from
%                              the eci frame to the body frame
%     w_b_ib                 - Angular velocity of the ECI frame with
%                              respect to the body frame, expressed in the
%                              body frame
%     B_body                 - Estimated magnetic field expressed in body frame
%     mtq_max                - Maximum dipole provided by each Magnetorquer
%     lim_dz                 - Absolute value of the limits of the deadzone [revolutions/min]
%     AngVel_rw_radps_cur    - RW angular velocity [rad/sec]
%     AngVel_rw_rpm_cur      - RW angular velocity [revolutions/min]
%     acceleration_rw_old    - RW acceleration
%     init_AngVel_dz         - Initial RW angular velocity when entering deadzone
%     init_accel_dz          - Initial RW acceleration when entering deadzone
%     timeflag_dz            - Counter which indicates the time present in the deadzone
%     rw_max_torque          - Maximum torque provided by each Magnetorquer
%     B_body_real            - Real magnetic field expressed in body frame
%     time                   - Current timestep
%     sun_vector_eci         - Unit vector from satellite to sun expressed in ECI frame
%     known_rm               - Estimated constant residual magnetic dipole
%     Reaction wheel deadzone behavior thresholds (in rpm/sec):
%     const1_accel           - Threshold for Case 1
%     const2_accel           - Threshold for Case 2
%     const3_accel           - Threshold for Case 3
%     const4_accel           - Threshold for Case 4
%     AngVel_rw_lim          - Angular velocity limit for RW desaturation
%     sun_desired            - Desired sun vector
%
%   Outputs:
%     torque                 - Total applied torque
%     T_rw                   - Torque provided by the RW
%     T_magnetic_effective   - Torque provided by the MTQs
%     V_rw                   - Voltage applied on RW
%     I_rw                   - Amperage applied on RW
%     P_thermal_rw           - Power consumed from RW
%     AngVel_rw_rpm_new      - Next RW angular velocity [revolutions/min]
%     AngVel_rw_radps_new    - Next RW angular velocity [rad/sec]
%     acceleration_rw_cur    - Current RW acceleration
%     rw_ang_momentum        - RW angular momentum
%     init_AngVel_dz         - Initial RW angular velocity when entering deadzone
%     init_accel_dz          - Initial RW acceleration when entering deadzone
%     V_mtq                  - Voltage applied on MTQs
%     I_mtq                  - Amperage applied on MTQs
%     P_thermal_mtq          - Power consumed from MTQs
%     timeflag_dz            - Counter which indicates the time present in the deadzone
%     q_sb                   - Quaternion that expresses the rotation from
%                              the sun-pointing frame to the body frame
% ======================================================================== %

function [torque, T_rw, T_magnetic_effective, V_rw, I_rw, P_thermal_rw, AngVel_rw_rpm_new, AngVel_rw_radps_new, ...
        acceleration_rw_cur, rw_ang_momentum, init_AngVel_dz, init_accel_dz, V_mtq, I_mtq, P_thermal_mtq, ...
        timeflag_dz, M, q_sb] = ...
        PD_Sun_Pointing(q_desired, q_eci_body, w_b_ib, B_body, mtq_max, ...
        lim_dz, AngVel_rw_radps_cur, AngVel_rw_rpm_cur, acceleration_rw_old, init_AngVel_dz, ...
        init_accel_dz, timeflag_dz, rw_max_torque, B_body_real, time, sun_vector_eci, known_rm,const1_accel,const2_accel,const3_accel,const4_accel,AngVel_rw_lim,sun_desired)

    global T_rw_data;
    global T_magnetic_data;
    global flag;
    global Jw;
    global Max_RW_torq;

    Kp_gain = 8e-03 * diag([1, 3, 1]);
    Kd_gain = 2e-01 * diag([1, 1, 1]);

    q_sb = q_sun_body(sun_vector_eci, q_eci_body,sun_desired);

    q_error = quatProd(quatconj(q_desired), q_sb);
    T_commanded = -sign(q_error(1)) * Kp_gain * q_error(2:4) - Kd_gain * w_b_ib;

    b_hat = B_body / norm(B_body);
    T_rw = [0; 0; 1] * (B_body' * T_commanded) / B_body(3);
    T_magnetic = skew(b_hat)' * skew(b_hat) * (T_commanded - T_rw);
    M = skew(B_body) * T_magnetic / (B_body' * B_body);

    %%  Saturation of the MTQs

    [T_magnetic, T_rw] = mtq_saturation(T_magnetic, T_rw, T_commanded, B_body, M, mtq_max, known_rm);

    M = -cross(T_magnetic, B_body) / (norm(B_body))^2;
    M = M - known_rm';
    T_magnetic_effective = cross(M, B_body_real);

    %%  Saturation of the RW

    if time > 1
        [T_magnetic_effective, T_rw] = ...
            rw_saturation(T_magnetic_effective, T_rw, acceleration_rw_old, AngVel_rw_rpm_cur, B_body,mtq_max,AngVel_rw_lim);

        if T_rw(3) > Max_RW_torq
            T_rw(3) = Max_RW_torq;
        elseif T_rw(3) < -Max_RW_torq
            T_rw(3) = -Max_RW_torq;
        end
    end

    %% Calculation of V_rw, I_rw, P_Rw in case of no-zero crossing

    if timeflag_dz == 0 && flag == 0
        [V_rw, I_rw, P_thermal_rw, AngVel_rw_radps_new, acceleration_rw_cur, T_rw_total] = ...
            rw_model(T_rw(3), AngVel_rw_radps_cur);
        AngVel_rw_rpm_new = 30 / pi * AngVel_rw_radps_new;
    end

    %%  Deadzone, are you here?
    if abs(AngVel_rw_rpm_cur) <= lim_dz && abs(T_rw(3)) > 9e-7
        if timeflag_dz == 0
            init_AngVel_dz = AngVel_rw_rpm_cur;
            init_accel_dz = acceleration_rw_cur;
        end
        [V_rw, I_rw, P_thermal_rw, AngVel_rw_rpm_new, acceleration_rw_cur, T_rw(3), timeflag_dz, init_accel_dz] = ...
            rw_deadzone(AngVel_rw_rpm_cur, timeflag_dz, init_accel_dz, init_AngVel_dz,const1_accel,const2_accel,const3_accel,const4_accel);
        AngVel_rw_radps_new = pi / 30 * AngVel_rw_rpm_new;

        T_magnetic = skew(b_hat)' * skew(b_hat) * (T_commanded - T_rw);
        M = skew(B_body) * T_magnetic / (B_body' * B_body);
        M = mtq_scaling(M, mtq_max);
        T_magnetic_effective = cross(M, B_body);
    else
        if timeflag_dz ~= 0
            timeflag_dz = 0;
            flag = 1;
        end
    end

    if timeflag_dz == 0 && flag == 1
        [V_rw, I_rw, P_thermal_rw, AngVel_rw_radps_new, acceleration_rw_cur] = ...
            rw_model(T_rw(3), AngVel_rw_radps_cur);
        AngVel_rw_rpm_new = 30 / pi * AngVel_rw_radps_new;
        flag = 0;
    end

    %% If the angular velocity is too small, it's actually 0

    if abs(AngVel_rw_rpm_new) < 0.5
        AngVel_rw_rpm_new = 0;
    end

    %%
    rw_ang_momentum = Jw * AngVel_rw_radps_new;

    %% Calculate V, I, P of MTQ's

    [V_mtq, I_mtq, P_thermal_mtq] = mtq_model(M);

    %% T_rw scaling

    if T_rw(3) > rw_max_torque
        T_rw(3) = rw_max_torque;
    elseif T_rw(3) < -rw_max_torque
        T_rw(3) = -rw_max_torque;
    end

    T_rw_data = [T_rw_data, T_rw];
    T_magnetic_data = [T_magnetic_data, T_magnetic_effective];
    torque = T_magnetic_effective + T_rw;

end