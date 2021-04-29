function [rotated_vector] = rotate_vector(quaternion, vector)
    temp = quatProd( quatconj(quaternion'), quatProd([0; vector], quaternion) );
    rotated_vector = temp(2:4);
end

