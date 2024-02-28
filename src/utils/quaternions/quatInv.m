%% Script for calculating the inverse of a quaternion
%
%  Input    :
%  Q        - quaternion
%
%  Output   :
%  invQ     - inverse of q
%
%
%  Q must have the form Q = [n e] where n is the scalar 
%  and e the vector part.
%  Note: Q does not need to be a unit quaternion.

function invQ = quatInv(Q)

invQ = zeros(size(Q));

n = (Q'*Q);
if (n < 1e-16)
    invQ = zeros(4,1);
    return;
end

invQ(1) = Q(1);
invQ(2:4) = -Q(2:4);

invQ = invQ / n;

end

