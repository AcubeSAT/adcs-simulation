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