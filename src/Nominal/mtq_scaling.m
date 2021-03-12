function [output] = mtq_scaling(M, mtq_max1, mtq_max2, mtq_max3)

output = M;
max1 = mtq_max1;  % Setting the maximum value that the MTQs can produce.
max2 = mtq_max2;
max3 = mtq_max3;

%% If M exceeds max value, then scale M to max
if(M(1) > max1)
    output(1) = max1;
end
if(M(2) > max2)
    output(2) = max2;
end
if(M(3) > max3)
    output(3) = max3;
end

%% If M exceeds -max value, then scale M to -max
if(M(1) < -max1)
    output(1) = -max1;
end
if(M(2) < -max2)
    output(2) = -max2;
end
if(M(3) < -max3)
    output(3) = -max3;
end

end