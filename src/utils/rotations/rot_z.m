%% Script for creating a 3x3 rotation matrix around z axis
%
%  Inputs   :
%  theta    - angle of rotation in degrees
%
%  Outputs  :
%  Rz       - rotation matrix

function Rz = rot_z(theta)

    Rz = [cosd(theta), -sind(theta), 0; sind(theta), cosd(theta), 0; 0, 0, 1];

end