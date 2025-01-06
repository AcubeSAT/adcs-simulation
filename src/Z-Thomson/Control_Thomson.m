
function [Bdot_body,torq,V_mtq, I_mtq, P_thermal_mtq,M] = ...
        Control_Thomson(w_b_ib, B2,B, mtq_max, Kd, Ks, w_ref)

Bx=acos(B(1)/norm(B));
Bx2=acos(B2(1)/norm(B2));

By=acos(B(2)/norm(B));
By2=acos(B2(2)/norm(B2));

Bz=acos(B(3)/norm(B));
Bz2=acos(B2(3)/norm(B2));

Bdot_body = ([Bx2;By2;Bz2] - [Bx;By;Bz]) / 0.1;
  

M(3,1)= Kd * Bdot_body(3);
if (abs(B(2))>abs(B(1)))
    M(1,1)= Ks*(abs(w_b_ib(3))-w_ref)*sign(B(2));
    M(2,1)= 0;
elseif (abs(B(2))<abs(B(1)))
    M(1,1)= 0;
    M(2,1)=-Ks*(abs(w_b_ib(3))-w_ref)*sign(B(1));
end  

M = mtq_scaling(M, mtq_max);
T_magnetic_effective = cross(M, B);
torq = T_magnetic_effective;

%% Calculate V, I, P of MTQ's
[V_mtq, I_mtq, P_thermal_mtq] = mtq_model(M);
       
end
   