%% ========================================================================
% Wahba's Problem
%% Inputs:
% sun_body = the sun vector in the body frame
% mgn_body = the magnetic field vector in the body frame
% sun_ref = the sun vector in the reference frame, in our design the ECI
% mgn_ref = the magnetic field vector in the reference frame, in our design
% the ECI frame.
%% NOTE
% The Wahba's Problem requires two pairs of vectors in two difference frames
% to calculate the quaternion from one frame to another. In our design
% these two frames are the Body of the Satellite and the ECI frame. However
% this function can be used for any other frames as well
%% Outputs:
% quat_r = the calculated quaternion from the reference frame (ECI) to the
% second frame (Body)
% Ropt_r = the Direction Cosine Matrix between these frames.
%%
% This function solves the Wahba's Probleb using the SVD Algorithm. The
% attitude matrix A is calculated using the pairs of vectors in the
% difference frames. In some applications, the vectors are weighted, i.e. 
% the calculation of A is based more on one of the two vectors, when the 
% second is not so trustworthy, but in our approach we consider that both 
% sensors are equally reliable, hence the weighs are 1.  The SVD algorithm 
% splits matrix A into three new matrices, two orthogonal (U,V) and one 
% diagonal. From these matrices we can find the Direction Cosine Matrix and 
% then calculate the normalised  quaternion. In our design, we use this 
% function only to find an initial state, which will be later used in the 
% Kalman Filter, since the performance of the Filter is significantly 
% better than the Wahba's Problem. 
%% 
% For deeper understanding of the Wahba's Problem, a term important in
% Attitude Determination, you are suggested to study Markley and Crassidis:
% Fundementals of Spacecraft Attitude Determination and Control.
% =========================================================================



function [quat_r, Ropt_r]=wahba(sun_body,mgn_body,sun_ref,mgn_ref)


mgn_ref_vec=mgn_ref'/norm(mgn_ref);
sun_ref_vec=sun_ref'/norm(sun_ref);

sun_wieght = 1;
mgn_weight=1;
A=zeros(3,3);

mgn_body=mgn_body'/norm(mgn_body);
sun_body=sun_body'/norm(sun_body);
for i=1:3
    A(i,1)=sun_wieght*sun_body(1)*sun_ref_vec(i)+mgn_weight*mgn_body(1)*mgn_ref_vec(i);
    A(i,2)=sun_wieght*sun_body(2)*sun_ref_vec(i)+mgn_weight*mgn_body(2)*mgn_ref_vec(i);
    A(i,3)=sun_wieght*sun_body(3)*sun_ref_vec(i)+mgn_weight*mgn_body(3)*mgn_ref_vec(i);
end

[U,S,V]=svd(A);

d=det(U)*det(V);

if d<0
    U(1,3)=-U(1,3);
    U(2,3)=-U(2,3);
    U(3,3)=-U(3,3);
end

Ropt_r=U*V';
quat_r = d2q(Ropt_r);                 
end