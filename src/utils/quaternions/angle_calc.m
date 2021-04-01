function [output] = angle_calc(x, y)
% This function calculates the cosine of two 1x3 vectors, x and y.

x_mag = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
y_mag = sqrt(y(1)^2 + y(2)^2 + y(3)^2);

XY_dot = x(1)*y(1) + x(2)*y(2) + x(3)*y(3);

output = XY_dot/(x_mag*y_mag);

end




