%% ======================================================================== %%
%   Function for scaling the desired magnetic dipole, in case it exceeds   
%   the maximum value that can be provided from each magnetorquer.
%   
%   Input: Desired Dipole, 
%          Maximum dipole of each magnetorquer (3 values)   
%   Ouput: Scaled Magnetic Dipole
% ========================================================================= %

function [scaled_dipole] = mtq_scaling(M, mtq_max1, mtq_max2, mtq_max3)

scaled_dipole = M;
max1 = mtq_max1; 
max2 = mtq_max2;
max3 = mtq_max3;

%% If M exceeds max value, then scale M to max
if(M(1) > max1)
    scaled_dipole(1) = max1;
end
if(M(2) > max2)
    scaled_dipole(2) = max2;
end
if(M(3) > max3)
    scaled_dipole(3) = max3;
end

%% If M exceeds -max value, then scale M to -max
if(M(1) < -max1)
    scaled_dipole(1) = -max1;
end
if(M(2) < -max2)
    scaled_dipole(2) = -max2;
end
if(M(3) < -max3)
    scaled_dipole(3) = -max3;
end

end