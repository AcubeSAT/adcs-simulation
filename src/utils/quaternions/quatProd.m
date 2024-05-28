%% Script for calculating the product of 2 quaternions
%
%  Inputs   :
%  Q1       - quaternion
%  Q2       - quaternion
%
%  Output   :
%  Q12      - product of Q1 and Q2
%
%
%  Each element of Q1 and Q2 must be a real number
%  Q1 and Q2 must have the form Q = [n e] where n
%  is the scalar and e the vector %   part.
%  Note: Quaternion multiplication is not commutative

function Q12 = quatProd(Q1, Q2)

    Q1 = Q1(:);
    Q2 = Q2(:);

    n1 = Q1(1);
    e1 = Q1(2:4);

    n2 = Q2(1);
    e2 = Q2(2:4);

    Q12 = zeros(4, 1);
    Q12(1) = n1 * n2 - e1' * e2;
    Q12(2:4) = n1 * e2 + n2 * e1 + cross(e1, e2);

end
