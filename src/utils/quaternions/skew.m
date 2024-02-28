%% Script for creating the skew symmetric form of a vector
%
%  Input    :
%  x        - 1x3 vector
%
%  Output   : 
%  X        - skew symmetric form of x

function X = skew(x)

X = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

end
