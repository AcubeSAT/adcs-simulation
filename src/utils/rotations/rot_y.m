%% Script for creating a 3x3 rotation matrix around y axis
%
%  Inputs   :
%  theta    - angle of rotation in degrees
%
%  Outputs  :
%  Ry       - rotation matrix


function Ry = rot_y(theta)

    Ry = [cosd(theta), 0, sind(theta); 0, 1, 0; -sind(theta), 0, cosd(theta)];

end