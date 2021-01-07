function [quaternion] = angles_to_quaternion(x,y,z)
    %x,y,z angles in each axis, in rads
    Rotation_x = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
    Rotation_y = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)]; %sign in sin according to Wikipedia
    Rotation_z = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];
    Rotation_final = Rotation_z*Rotation_y*Rotation_x;
    quaternion = d2q(Rotation_final);
end

