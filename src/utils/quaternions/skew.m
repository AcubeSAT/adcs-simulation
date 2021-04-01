function X = skew(x)
%% Input:
% x = vector 1x3
%% Output 
% X = vector, 1x3, the skew symmetric form of x

    X = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
end
