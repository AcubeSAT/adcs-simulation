function v_rot = quatLog(Q)
%   v_rot = quatLog(Q) calculates the logarithmic map of quaternion Q,
%   which is the rotational velocity, v_rot, that in unit time causes
%   rotation equal to Q.

    n = Q(1);
    e = Q(2:4);
    norm_e = norm(e);
    
    if (norm_e > 1e-16)
        % use real(acos(η)) because for n close but greater than 1, acos(n) becomes complex number
%         theta = 2*real(acos(n));
        theta = 2*real(atan2(norm_e,n));
        v_rot = theta*e/norm_e;
    else
        v_rot = zeros(size(e));
    end

end



