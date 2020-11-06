
function [quat_r Ropt_r]=wahba(sun,mgn,sun_ref,mgn_ref)

%% Create reference, observation and weight vectors


mgn_ref_vec=mgn_ref'/norm(mgn_ref);
sun_ref_vec=sun_ref'/norm(sun_ref);

sun_gain = 1;
mgn_gain=1;
g=0;
ga=1-sun_gain;
gm=1-mgn_gain;
A=zeros(3,3);

mgn=mgn'/norm(mgn);
sun=sun'/norm(sun);
for i=1:3
    A(i,1)=sun_gain*sun(1)*sun_ref_vec(i)+mgn_gain*mgn(1)*mgn_ref_vec(i);
    A(i,2)=sun_gain*sun(2)*sun_ref_vec(i)+mgn_gain*mgn(2)*mgn_ref_vec(i);
    A(i,3)=sun_gain*sun(3)*sun_ref_vec(i)+mgn_gain*mgn(3)*mgn_ref_vec(i);
end

[U,S,V]=svd(A);

d=det(U)*det(V);
S(2,2)+d*S(3,3);

if d<0
    U(1,3)=-U(1,3);
    U(2,3)=-U(2,3);
    U(3,3)=-U(3,3);
end

Ropt_r=U*V';
quat_r = d2q(Ropt_r);                 
end