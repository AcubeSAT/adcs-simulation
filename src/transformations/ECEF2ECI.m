function [vec_eci, v_vec_eci]= ECEF2ECI(vec_ecef,v_vec_ecef,gst,vel)
%Now rotate the coordinates
CGAST = cos(-gst); SGAST = sin(-gst);
R=[CGAST SGAST 0; -SGAST CGAST 0; 0 0 1];
% temp=R(2,:);
% R(2,:)=R(1,:);
% R(1,:)=temp;
vec_eci=R*vec_ecef;
if vel == 1 
    %Apply rotation to convert velocity vector from ECI to ECEF coordinates
    OMEGAE = 7.29211586D-5;  %Earth rotation rate in rad/s
    temp=(CGAST^(2)-SGAST)
    v_vec_eci(1)=(v_vec_ecef(1)*CGAST-v_vec_ecef(2)+vel*OMEGAE*vec_eci(1)-CGAST*vel*OMEGAE*vec_eci(2))/temp;
    v_vec_eci(2)=(v_vec_ecef(1)*SGAST*CGAST+v_vec_ecef(2)*(temp-SGAST)+vel*OMEGAE*vec_eci(1)*(temp-1)-vel*OMEGAE*vec_eci(2)*SGAST*CGAST)/(CGAST*(temp));
    v_vec_eci(3)= v_vec_ecef(3);
else
    v_vec_eci=zeros(3,1);
end
v_vec_eci=R';
end