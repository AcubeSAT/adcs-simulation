% -----------------------------------------------------------------------------
%
%                              procedure angles_to_quaternion
%
%  this procedure creates a quaternion from the euler angles
%
%  Inputs         :
%   x             - Euler Angles in X axis
%   y             - Euler Angles in Y axis
%   z             - Euler Angles in Z axis
%
%  Outputs        :
%   quaternion    - quaternion from the Euler angles
%  ----------------------------------------------------------------------------*/


function [quaternion] = angles_to_quaternion(x,y,z)

% Construct rotation matrix
Rotation_x = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
Rotation_y = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)];
Rotation_z = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];
Rotation_final = Rotation_z*Rotation_y*Rotation_x;

% Rotation matrix to quaternion conversion
quaternion = d2q(Rotation_final);

end

