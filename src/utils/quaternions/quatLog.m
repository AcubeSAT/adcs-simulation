%% Script that calculates the logarithmic map of a quaternion
% 
%  Input    :
%  Q        - quaternion
%
%  Output   :
%  v_rot    - rotational velocity that in unit time causes rotation equal to Q
%
%
%  The if statement is used, because for small values of the norm_e,
%  the function atan is not well defined.

function v_rot = quatLog(Q)

n = Q(1);
e = Q(2:4);
norm_e = norm(e);
    
if (norm_e > 1e-16)
    theta = 2*real(atan2(norm_e,n));
    v_rot = theta*e/norm_e;
else
    v_rot = zeros(size(e));
end

end



