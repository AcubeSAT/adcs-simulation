function Q = e2q_ypr(euler)
    % Converts Euler Angles to representative Quaternions
    %
    % Copyright (C) 2013 - Pieter V. Reyneke, South Africa
    %
    % Open Source licence, use whereever you like but at own risk, keep this
    % complete copyright notice in the code and please send me updates and
    % report errors to pieter.reyneke@gmail.com. Will give recognition, if
    % requested, after any valid comment leading to an update was recieved.
    %
    % Author: PV Reyneke

    if (nargin < 1)
        n = 10000;
        by = 1 .* (rand(n, 1) .* 2 * pi - pi);
        bp = 1 .* (asin(rand(n, 1).*2-1));
        br = 1 .* (rand(n, 1) .* 2 * pi - pi);
        eo = [by, bp, br];
        bq = e2q_ypr(eo);
        eq = q2e_ypr(bq);
        display(num2str(sum(ssqr(eo-eq, 2))));
        return
    end

    spsi_2 = sin(euler(:, 1)/2);
    sthe_2 = sin(euler(:, 2)/2);
    sphi_2 = sin(euler(:, 3)/2);
    cpsi_2 = cos(euler(:, 1)/2);
    cthe_2 = cos(euler(:, 2)/2);
    cphi_2 = cos(euler(:, 3)/2);

    Q(:, 1) = (cphi_2 .* cthe_2 .* cpsi_2) + (sphi_2 .* sthe_2 .* spsi_2);
    Q(:, 2) = (sphi_2 .* cthe_2 .* cpsi_2) - (cphi_2 .* sthe_2 .* spsi_2);
    Q(:, 3) = (cphi_2 .* sthe_2 .* cpsi_2) + (sphi_2 .* cthe_2 .* spsi_2);
    Q(:, 4) = (cphi_2 .* cthe_2 .* spsi_2) - (sphi_2 .* sthe_2 .* cpsi_2);

    % Exersize constraint
    dq = 1 - sum(Q.^2, 2);
    Q = ((0.5 * repmat(dq, [1, 4])) + 1) .* Q;

end