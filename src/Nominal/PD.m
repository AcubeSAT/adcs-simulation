function  [torque, T_rw, T_magnetic_effective, V_rw, I_rw, P_thermal_rw, AngVel_rw_rpm_new, AngVel_rw_radps_new,...
            acceleration_rw_cur, rw_ang_momentum, init_AngVel_dz, init_accel_dz, V_mtq, I_mtq, P_thermal_mtq, ...
                timeflag_dz] = ...
                    PD(Kp_gain, Kd_gain, q_desired , q_orbit_body , w_b_ib , B_body , eclipse, mtq_max, ...
                        lim_dz, AngVel_rw_radps_cur, AngVel_rw_rpm_cur, acceleration_rw_old, init_AngVel_dz, ...
                            init_accel_dz,timeflag_dz,rw_max_torque,B_body_real,time)
    
    global T_rw_data;
    global T_magnetic_data;
    global flag;
    global Jw;
    global Max_RW_torq;
    
%     Kp_gain= 1e-04*diag([5 10 5]);                                                  % Calculate gain (wrt Markley-Crassidis)
%     Kd_gain= 1e-02*diag([1 1 1]);
    w_o_io = [0;0.00110808802079241;0];
    %w_o = 0.00113308802079241;
    
%     if q_orbit_body(1)<0
%         q_orbit_body = -q_orbit_body;
%     end
    
  %  R_OB = quat2dcm(q_orbit_body'); % Calculating the transformation matrix from orbit to body frame

    q_w_b_io = quatProd(quatconj(q_orbit_body') ,quatProd([0;w_o_io],q_orbit_body));
    w_b_io = q_w_b_io(2:4);

    w_b_ob = w_b_ib - w_b_io; % Calculating angular rate of satellite relative to ECI frame

    q_error=quatProd(conj(q_desired),q_orbit_body);
    T_commanded = -sign(q_error(1))*Kp_gain*q_error(2:4) - Kd_gain*w_b_ob;

    b_hat=B_body/norm(B_body); 
    T_rw =[0;0;1]*(B_body'*T_commanded)/B_body(3);
    T_magnetic = skew(b_hat)'*skew(b_hat) * (T_commanded-T_rw);   
    M=skew(B_body)*T_magnetic/(B_body'*B_body);

    % Calculating gains in case of saturation

    M2 = (1/norm(T_magnetic))*M;
    Tm = T_magnetic/norm(T_magnetic);
    Tw = T_rw/norm(T_rw);
    Kma = (Tm'-(Tw'*Tm)*Tw')*T_commanded/(1-(Tw'*Tm)^2);
    Kwa = (Tw'-(Tm'*Tw)*Tm')*T_commanded/(1-(Tm'*Tw)^2);
    if max(abs(M)) > mtq_max
    Kma_s = min(abs(mtq_max./M2));
    Kwa_s = Kma_s*Kwa/Kma;
    Ms = Kma_s*M2;
    T_magnetic = skew(B_body)'*Ms;
    T_rw = Kwa_s.*Tw;
    else
    T_magnetic = Kma.*Tm;
    T_rw = Kwa.*Tw;
    end

    M = -cross(T_magnetic,B_body)/(norm(B_body))^2;
    %   M = mtq_scaling(M, mtq_max);
    %   M=M-rm;
    T_magnetic_effective = cross(M,B_body);

    %%  Saturation of the RW
   
    if time > 1 
    [T_magnetic_effective, T_rw] = ...
        rw_saturation(T_magnetic_effective, T_rw, acceleration_rw_old, AngVel_rw_rpm_cur, B_body);
    
    if T_rw(3) > Max_RW_torq
        T_rw(3) = Max_RW_torq;
    elseif T_rw(3) < -Max_RW_torq
        T_rw(3) = -Max_RW_torq;
    end
    
    T_commanded = T_magnetic_effective + T_rw;
   end

    %% Calculation of V_rw, I_rw, P_Rw in case of no-zero crossing

    if timeflag_dz == 0 && flag == 0
        [V_rw, I_rw, P_thermal_rw, AngVel_rw_radps_new, acceleration_rw_cur, T_rw_total] = ...
                                                rw_model(T_rw(3), AngVel_rw_radps_cur);
        AngVel_rw_rpm_new = 30/pi * AngVel_rw_radps_new;  
    end

    %%  Deadzone, are you here?
    if abs(AngVel_rw_rpm_cur) <= lim_dz && abs(T_rw(3)) > 9e-7
        if timeflag_dz == 0
            init_AngVel_dz = AngVel_rw_rpm_cur;
            init_accel_dz = acceleration_rw_cur;
        end
        [V_rw, I_rw, P_thermal_rw, AngVel_rw_rpm_new, acceleration_rw_cur, T_rw(3), timeflag_dz, init_accel_dz] = ...
            rw_deadzone(AngVel_rw_rpm_cur, timeflag_dz, init_accel_dz, init_AngVel_dz);
        AngVel_rw_radps_new = pi/30 * AngVel_rw_rpm_new;

        T_magnetic = skew(b_hat)'*skew(b_hat) * (T_commanded - T_rw);   
        M = skew(B_body)*T_magnetic/(B_body'*B_body);
        M = mtq_scaling(M, mtq_max);
        T_magnetic_effective = cross(M,B_body);
    else
       if timeflag_dz ~= 0
           timeflag_dz = 0;
           flag = 1;
       end
    end 

    if timeflag_dz == 0 && flag == 1
        [V_rw, I_rw, P_thermal_rw, AngVel_rw_radps_new, acceleration_rw_cur] = ...
                                                rw_model(T_rw(3), AngVel_rw_radps_cur);
        AngVel_rw_rpm_new = 30/pi * AngVel_rw_radps_new;
        flag = 0;
    end

    %% If the angular velocity is too small, it's actually 0

    if abs(AngVel_rw_rpm_new) < 0.5
            AngVel_rw_rpm_new = 0;
    end

    %%
    rw_ang_momentum = Jw * AngVel_rw_radps_new;

    %%  Calculate V, I, P of MTQ's

    [V_mtq, I_mtq, P_thermal_mtq] = mtq_model(M);

    %% T_rw scaling 
    if T_rw(3) > rw_max_torque
        T_rw(3) = rw_max_torque;
    elseif T_rw(3) < -rw_max_torque
        T_rw(3) = -rw_max_torque;
    end
    
    T_rw_data = [T_rw_data T_rw];
    T_magnetic_data= [T_magnetic_data T_magnetic_effective];
    torque = T_magnetic_effective + T_rw;
   
end