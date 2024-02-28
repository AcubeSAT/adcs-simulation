%% Script that calculates the exponential map of the rotational velocity 
%
%  Input     :
%  v_rot     - rotational velocity
%
%  Output    :
%  Q         - quaternion that causes rotation equal to that of v_rot in unit time
%
%  Algorithm :
%  v_rot = k*theta, where k is the normalized angle of rotation
%  norm(v_rot) = |theta|
%  v_rot/norm(v_rot) = k * sign(sin(theta))
%
%
%  The if statement is used to avoid division with small values of norm_v_rot

function Q = quatExp(v_rot)
  
norm_v_rot = norm(v_rot);
theta = norm_v_rot;

Q = zeros(4,1);

if (norm_v_rot > 1e-16)
    Q(1) = cos(theta/2);
    Q(2:4) = sin(theta/2) * v_rot / norm_v_rot;
else
    Q = [1 0 0 0]';
end

end


