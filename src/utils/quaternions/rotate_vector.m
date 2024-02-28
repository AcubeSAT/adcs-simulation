%% Script for rotating a vector by a quaternion
%
%  Inputs           :
%  quaternion       - quaternion representing the rotation
%  vector           - 3x1 vector to rotate
%
%  Output           :
%  rotated_vector   - rotated vector
%
%
%  Q must be a unit quaternion.

function [rotated_vector] = rotate_vector(quaternion, vector)

% quaternion_norm = sqrt(quaternion' * quaternion);
% quaternion = quaternion / quaternion_norm;

temp = quatProd( quatconj(quaternion'), quatProd([0; vector], quaternion) );

rotated_vector = temp(2:4);

end

