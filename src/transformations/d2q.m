
function q=d2q(dcm)
%D2Q Direction Cosine Matrix to normalized Quaternion
%	q=d2q(dcm), with:
%
%	dcm	Direction cosine matrix
%	q	Quaternion q=[q0,q1,q2,q3]'
%
% Copyright (C) 2013 - Pieter V. Reyneke, South Africa
%
% Open Source licence, use whereever you like but at own risk, keep this
% complete copyright notice in the code and please send me updates and 
% report errors to pieter.reyneke@gmail.com. Will give recognition, if
% requested, after any valid comment leading to an update was recieved. 
%
% Author: PV Reyneke
if (nargin<1)
    n = 1;
    y = rand(n,1).*2*pi-pi;
    p = asin(rand(n,1).*2-1);
    r = rand(n,1).*2*pi-pi;
    qqq = e2q_ypr([y p r]);
    dcm = q2d(qqq)';
end
tr=dcm(1,1)+dcm(2,2)+dcm(3,3);
Pa=1+tr;
Pb=1+2*dcm(1,1)-tr;
Pc=1+2*dcm(2,2)-tr;
Pd=1+2*dcm(3,3)-tr;
q=zeros(1,4);
if (Pa>=Pb && Pa>=Pc && Pa>=Pd)
    q(1,1)=0.5*sqrt(Pa);
    q(1,2)=(dcm(3,2)-dcm(2,3))/4/q(1);
    q(1,3)=(dcm(1,3)-dcm(3,1))/4/q(1);
    q(1,4)=(dcm(2,1)-dcm(1,2))/4/q(1);
elseif (Pb>=Pc && Pb>=Pd)
    q(1,2)=0.5*sqrt(Pb);
    q(1,3)=(dcm(2,1)+dcm(1,2))/4/q(2);
    q(1,4)=(dcm(1,3)+dcm(3,1))/4/q(2);
    q(1,1)=(dcm(3,2)-dcm(2,3))/4/q(2);
elseif (Pc>=Pd)
    q(1,3)=0.5*sqrt(Pc);
    q(1,4)=(dcm(3,2)+dcm(2,3))/4/q(3);
    q(1,1)=(dcm(1,3)-dcm(3,1))/4/q(3);
    q(1,2)=(dcm(2,1)+dcm(1,2))/4/q(3);
else
    q(1,4)=0.5*sqrt(Pd);
    q(1,1)=(dcm(2,1)-dcm(1,2))/4/q(4);
    q(1,2)=(dcm(1,3)+dcm(3,1))/4/q(4);
    q(1,3)=(dcm(3,2)+dcm(2,3))/4/q(4);
end
if (q(1)<=0)
    q=-q;
end
if (nargin<1)
    [dcm q2d(q) (dcm-q2d(q)) dcm-e2d_ypr([y p r])]
end
