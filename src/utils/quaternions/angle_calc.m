% *************************************************************************
% @file           : angle_calc.m
% @brief          : Function that calculates the cosine of two 1x3 vectors
% *************************************************************************
% @attention
% 
% <h2><center>&copy; Copyright (c) Aristotle Space & Aeronautics Team.
% All rights reserved.</center></h2>
%  
% This project is licensed under the GNU General Public License v3.0.
%  
% *************************************************************************

function [output] = angle_calc(x, y)

x_mag = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
y_mag = sqrt(y(1)^2 + y(2)^2 + y(3)^2);

XY_dot = x(1)*y(1) + x(2)*y(2) + x(3)*y(3);

output = XY_dot/(x_mag*y_mag);

end




