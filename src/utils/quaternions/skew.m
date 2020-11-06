% *************************************************************************
% @file           : skew.m
% @brief          : Function that calculates the skew symmetric form of a
%                   1x3 vector.
% *************************************************************************
% @attention
% 
% <h2><center>&copy; Copyright (c) Aristotle Space & Aeronautics Team.
% All rights reserved.</center></h2>
%  
% This project is licensed under the GNU General Public License v3.0.
%  
% *************************************************************************

function X = skew(x)
    %Returns the skew symmetric form of x (1x3)
    X = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
end
