%% Script for creating a 3x3 rotation matrix around x axis
%
%  Inputs   :
%  theta    - angle of rotation in degrees
%
%  Outputs  :
%  Rx       - rotation matrix


function Rx = rot_x(theta)

    Rx = [1, 0, 0; 0, cosd(theta), -sind(theta); 0, sind(theta), cosd(theta)];

end