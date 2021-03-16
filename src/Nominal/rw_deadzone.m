%% ========================================================================%%
%   Function dedicated for the behaviour of the RW during deadzone.
%   Modelled using data from the component's datasheet. It will be updated
%       after specified testing for the RW.
%   
%   Inputs:  Current angular velocity of the RW, 
%            Counter which indicates the time present in the deadzone,
%            Initial RW acceleration when entering deadzone,
%            Initial RW angular velocity when entering deadzone
%   Outputs: Voltage applied on RW, 
%            Amperage applied on RW, 
%            Power consumed from RW,  
%            Next angular velocity of the RW, 
%            Acceleration of the RW,
%            Next RW Torque, 
%            Updated counter, 
%            Initial RW acceleration when entering deadzone
% =========================================================================%

function [V_rw, I_rw, P_thermal_rw, AngVel_rw_new, accel_rw, T_rw_new, timeflag_dz, init_accel_dz] = ...
                rw_deadzone(AngVel_rw, timeflag_dz, init_accel_dz, initAngVel_dz)

%   The changes of the acceleration above happen inside a timeperiod of 
%   500 * 2 / 100 = 11sec if abs(accel) == 25.
%   400 * 2 / 100 = 8sec if abs(accel) == 50.
%   350 * 2 / 100 = 7sec if abs(accel) == 100.
%   300 * 2 / 100 = 7sec if abs(accel) == 200.
    
    global Km;
    global Rb;
    global Kv;
    global Jw;
    global c_friction;
    global b_friction;
    global Max_RW_torq;

    const1_accel = 25;
    const2_accel = 50;
    const3_accel = 100;
    const4_accel = 200;
    case1 = 0;
    case2 = 0;
    case3 = 0;
    case4 = 0;
    timestep = 0.1; 
    accel_rw = init_accel_dz;
    
    
    if timeflag_dz == 0
        if abs(init_accel_dz) > 0 && abs(init_accel_dz) < 37
            case1 = 1;
            divider = const1_accel/init_accel_dz;
        elseif abs(init_accel_dz) >= 37 && abs(init_accel_dz) < 75
            case2 = 1;
            divider = const2_accel/init_accel_dz;
        elseif abs(init_accel_dz) >= 75 && abs(init_accel_dz) < 150
            case3 = 1;
            divider = const3_accel/init_accel_dz;
        else
            case4 = 1;
            divider = const4_accel/init_accel_dz;
        end
    end
    
    timeflag_dz = timeflag_dz + timestep;    % The flag_rw indicates the time that has passed since the
                                             % start of the zero-crossing
        
    %% The case of abs(accel_rw) ~= 25 rpm/sec
    
    if case1 == 1
    
       check1_accel = 25/divider;
       check2_accel = 0;
       check3_accel = 300/divider;
       check4_accel = 250/divider;
       check5_accel = 75/divider;
       check6_accel = 75/divider;
       check7_accel = 25/divider;
       
       check1_time = 4.4;
       check2_time = 4.9;
       check3_time = 5.9;
       check4_time = 7.3;
       check5_time = 7.9;
       check6_time = 8.9;
       check7_time = 9.9;
       
       if abs(initAngVel_dz) < 10 && timeflag_dz == 0.1  % If the starting value of velocity is almost 0
           timeflag_dz = 0.1;
       end
       
       if timeflag_dz > 0 && timeflag_dz < 4.5
                accel_rw = sign(init_accel_dz)*check1_accel;
       elseif timeflag_dz >= 4.5 && timeflag_dz < 5
                accel_rw = sign(init_accel_dz)*(check1_accel - 50/divider * abs(timeflag_dz - check1_time));
       elseif timeflag_dz >= 5 && timeflag_dz < 6
                accel_rw = sign(init_accel_dz)*(check2_accel + 300/divider * abs(timeflag_dz - check2_time));
       elseif timeflag_dz >= 6 && timeflag_dz < 7.4
                accel_rw = sign(init_accel_dz)*(check3_accel - 35/divider * abs(timeflag_dz - check3_time));
       elseif timeflag_dz >= 7.4 && timeflag_dz < 8
                accel_rw = sign(init_accel_dz)*(check4_accel - 290/divider * abs(timeflag_dz - check4_time));
       elseif timeflag_dz >= 8 && timeflag_dz < 9
                accel_rw = sign(init_accel_dz)*(check5_accel + 0 * abs(timeflag_dz - check5_time));
       elseif timeflag_dz >= 9 && timeflag_dz < 10
                accel_rw = sign(init_accel_dz)*(check6_accel - 50/divider * abs(timeflag_dz - check6_time));
       elseif timeflag_dz >= 10 && timeflag_dz < 11
                accel_rw = sign(init_accel_dz)*(check7_accel + 0 * abs(timeflag_dz - check7_time));
       end
           
       if timeflag_dz == 11
          timeflag_dz = 0;
%           fprintf('The angular velocity after end of the zero-crossing is: %f', AngVel_rw);
       end
    end 
    
    %% The case of abs(accel_rw) ~= 50 rpm/sec
    
    if case2 == 1
    
       check1_accel = 50/divider;
       check2_accel = 0;
       check3_accel = 350/divider;
       check4_accel = 350/divider;
       check5_accel = 150/divider;
       check6_accel = 50/divider;
       
       check1_time = 2.3;
       check2_time = 2.7;
       check3_time = 3.9;
       check4_time = 5.1;
       check5_time = 5.5;
       check6_time = 6.9;
       
       if abs(initAngVel_dz) < 10 && timeflag_dz == 0.1  % If the starting value of velocity is almost 0
           timeflag_dz = 0.1;
       end
        
       if timeflag_dz > 0 && timeflag_dz < 2.4
            accel_rw = sign(init_accel_dz)*check1_accel;
       elseif timeflag_dz >= 2.4 && timeflag_dz < 2.8
            accel_rw = sign(init_accel_dz)*(check1_accel - 125/divider * abs(timeflag_dz - check1_time));
       elseif timeflag_dz >= 2.8 && timeflag_dz < 4
            accel_rw = sign(init_accel_dz)*(check2_accel + 275/divider * abs(timeflag_dz - check2_time));
       elseif timeflag_dz >= 4 && timeflag_dz < 5.2
            accel_rw = sign(init_accel_dz)*(check3_accel + 0 * abs(timeflag_dz - check3_time));
       elseif timeflag_dz >= 5.2 && timeflag_dz < 5.6
            accel_rw = sign(init_accel_dz)*(check4_accel - 500/divider * abs(timeflag_dz - check4_time));
       elseif timeflag_dz >= 5.6 && timeflag_dz < 7
            accel_rw = sign(init_accel_dz)*(check5_accel - 70/divider * abs(timeflag_dz - check5_time));
       elseif timeflag_dz >= 7 && timeflag_dz < 8
            accel_rw = sign(init_accel_dz)*(check6_accel + 0 * abs(timeflag_dz - check6_time));
       end
       
       if timeflag_dz == 8
          timeflag_dz = 0;
%           fprintf('The angular velocity after end of the zero-crossing is: %f', AngVel_rw);
       end
    end 
    
    %% The case of abs(accel_rw) ~= 100 rpm/sec
    
    if case3 == 1
    
       check1_accel = 80/divider;
       check2_accel = 400/divider;
       check3_accel = 450/divider;
       check4_accel = 250/divider;
       check5_accel = 250/divider;
       check6_accel = 100/divider;
       
       check1_time = 1.5;
       check2_time = 2.9;
       check3_time = 3.5;
       check4_time = 4.4;
       check5_time = 4.7;
       check6_time = 5.9;
       
       if abs(initAngVel_dz) < 10 && timeflag_dz == 0.1  % If the starting value of velocity is almost 0
           timeflag_dz = 0.1;
       end
       
       if timeflag_dz > 0 && timeflag_dz < 1.6
            accel_rw = sign(init_accel_dz)*check1_accel;
       elseif timeflag_dz >= 1.6 && timeflag_dz < 3
            accel_rw = sign(init_accel_dz)*(check1_accel + 230/divider * abs(timeflag_dz - check1_time));
       elseif timeflag_dz >= 3 && timeflag_dz < 3.6
            accel_rw = sign(init_accel_dz)*(check2_accel + 85/divider * abs(timeflag_dz - check2_time));
       elseif timeflag_dz >= 3.6 && timeflag_dz < 4.5
            accel_rw = sign(init_accel_dz)*(check3_accel - 220/divider * abs(timeflag_dz - check3_time));
       elseif timeflag_dz >= 4.5 && timeflag_dz < 4.8
            accel_rw = sign(init_accel_dz)*(check4_accel + 0 * abs(timeflag_dz - check4_time));
       elseif timeflag_dz >= 4.8 && timeflag_dz < 6
            accel_rw = sign(init_accel_dz)*(check5_accel - 125/divider * abs(timeflag_dz - check5_time));
       elseif timeflag_dz > 6 && timeflag_dz < 7
            accel_rw = sign(init_accel_dz)*(check6_accel + 0 * abs(timeflag_dz - check6_time));
       end
       
       if timeflag_dz == 7
          timeflag_dz = 0;
%           fprintf('The angular velocity after end of the zero-crossing is: %f', AngVel_rw);
       end
    end 
    
    %% The case of abs(accel_rw) ~= 200 rpm/sec
    
    if case4 == 1
    
       check1_accel = 180/divider;
       check2_accel = 50/divider;
       check3_accel = 600/divider;
       check4_accel = 400/divider;
       check5_accel = 200/divider;
       
       check1_time = 1.4 - 0.1;
       check2_time = 3 - 0.1;
       check3_time = 3.5 - 0.1;
       check4_time = 6 - 0.1;
       
       if abs(initAngVel_dz) < 10 && timeflag_dz == 0.1  % If the starting value of velocity is almost 0
           timeflag_dz = 0.1;
       end
       
       if timeflag_dz > 0 && timeflag_dz < 1.4
            accel_rw =sign(init_accel_dz)*(check1_accel - 95/divider * timeflag_dz);
       elseif timeflag_dz >= 1.4 && timeflag_dz < 3
            accel_rw = sign(init_accel_dz)*(check2_accel + 345/divider * abs(timeflag_dz - check1_time));
       elseif timeflag_dz >= 3 && timeflag_dz < 3.5
            accel_rw = sign(init_accel_dz)*(check3_accel - 400/divider * abs(timeflag_dz - check2_time));
       elseif timeflag_dz >= 3.5 && timeflag_dz < 6
            accel_rw = sign(init_accel_dz)*(check4_accel - 80/divider * abs(timeflag_dz - check3_time));
       elseif timeflag_dz >= 6 && timeflag_dz < 7
            accel_rw = sign(init_accel_dz)*(check5_accel + 0 * abs(timeflag_dz - check4_time));
       end
       
       if timeflag_dz == 7
          timeflag_dz = 0;
%           fprintf('The angular velocity after end of the zero-crossing is: %f', AngVel_rw);
       end
    end 
    
    
    %% Calculate the new T_rw and the V_rw, I_rw, P_rw

    T_rw_new = (pi/30)*accel_rw * Jw + b_friction * (pi/30) * AngVel_rw + c_friction * sign(AngVel_rw); %In Nm
    I_rw = T_rw_new / Km + Kv * (pi/30) * AngVel_rw / Rb;
    AngVel_rw_new = AngVel_rw + accel_rw * timestep;    %In rpm = rpm + (rpm/sec)*sec
    V_rw = I_rw * Rb;
    P_thermal_rw = I_rw * V_rw;
    
    if T_rw_new > Max_RW_torq
        T_rw_new = Max_RW_torq;
    elseif T_rw_new < -Max_RW_torq
        T_rw_new = -Max_RW_torq;
    end

end