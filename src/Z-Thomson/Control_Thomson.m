% ===============================================================================================
%  Implementation of the controller that calculates the desired torque
%  to be applied in order to achieve Z-Thomson
%
%  Inputs:
%   w_b_ib         -Angular velocity of the ECI frame with
%                   respect to the body frame, expressed in the
%                   body frame
%  B_body          -Estimated magnetic field expressed in body frame
%  B_body_2        -Second magnetometer measurement 
%  B_body_real     -Real magnetic field expressed in body frame
%  mtq_max         -Maximum dipole provided by each Magnetorquer
%  Kd              -Gain for Detumbling
%  Ks              -Gain for Thomson
%  w_ref           -Desired angular velocity ,around z-axis of body frame with
%                   respect to the ECI frame , expressed in the body frame 
%
%  Outputs:
%  Bdot_body       -Magnetic field Deviation 
%  torq            -Applied torque
%  V_mtq           -Voltage applied on mtq
%  I_mtq           -Amperage applied on mtq
%  P_thermal_mtq   -Power consumed from mtq 
%  M               -Magnetic dipole
% ===============================================================================================

function [Bdot_body,torq,V_mtq, I_mtq, P_thermal_mtq,M] = ...
        Control_Thomson(w_b_ib, B_body_2,B_body,B_body_real ,mtq_max, Kd, Ks, w_ref)

    Bx=acos(B_body(1)/norm(B_body));
    Bx2=acos(B_body_2(1)/norm(B_body_2));

    By=acos(B_body(2)/norm(B_body));
    By2=acos(B_body_2(2)/norm(B_body_2));

    Bz=acos(B_body(3)/norm(B_body));
    Bz2=acos(B_body_2(3)/norm(B_body_2));

    Bdot_body = ([Bx2;By2;Bz2] - [Bx;By;Bz]) / 0.1;
    

    M(3,1)= Kd * Bdot_body(3);
    if (abs(B_body(2))>abs(B_body(1)))
        M(1,1)= Ks*(abs(w_b_ib(3))-w_ref)*sign(B_body(2));
        M(2,1)= 0;
    elseif (abs(B_body(2))<abs(B_body(1)))
        M(1,1)= 0;
        M(2,1)=-Ks*(abs(w_b_ib(3))-w_ref)*sign(B_body(1));
    end  

    M = mtq_scaling(M, mtq_max);
    T_magnetic_effective = cross(M, B_body_real);
    torq = T_magnetic_effective;

    %% Calculate V, I, P of MTQ's
    [V_mtq, I_mtq, P_thermal_mtq] = mtq_model(M);
        
end
   