function [T_mtq_effective, T_rw] = rw_saturation(T_magnetic_effective, T_rw, accel_rw, AngVel_rw, B_body, mtq_max1, mtq_max2, mtq_max3)

    % The angular velocity of the Reaction Wheel is given by a sensor
    % placed on the RW.
    
    global A;
    global Jw;
    
    AngVel_rw_lim = 10000;
    T_mtq_effective = T_magnetic_effective;

    if abs(AngVel_rw) > AngVel_rw_lim && abs(T_rw(3))> 0
        T_added = [0; 0; A * Jw * accel_rw];
        T_magnetic = T_magnetic_effective + T_added;
        M = -cross(T_magnetic,B_body)/(norm(B_body))^2;
        M = mtq_scaling(M, mtq_max1, mtq_max2, mtq_max3);
        T_mtq_effective = cross(M,B_body);
        
        T_rw = T_rw - (T_mtq_effective - T_magnetic_effective);
    end
end

