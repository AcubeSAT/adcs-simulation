% *************************************************************************
% @file           : mtq_scaling.m
% @brief          : Function that scales the magnitude of the magnetorquers
%                   whenever they exceed the maximum value.
% *************************************************************************
% @attention
% 
% <h2><center>&copy; Copyright (c) Aristotle Space & Aeronautics Team.
% All rights reserved.</center></h2>
%  
% This project is licensed under the GNU General Public License v3.0.
%  
% *************************************************************************

function [output] = mtq_scaling(M, mtq_max)

output = M;
max = mtq_max;  % Setting the maximum value that the MTQs can produce.

%% If M exceeds max value, then scale M to max
for i = 1:3
    if(M(i) > max)
        output(i) = max;
    end
end

%% If M exceeds -max value, then scale M to -max
for i = 1:3
    if(M(i) < -max)
        output(i) = -max;
    end
end

end