function euler = q2e_ypr(quaternion)
    % Converts Quaternions to representative Euler Angles
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
        euler = [1, 0.5, 0.1];
        quaternion = e2q_ypr([-1, 0.5, -0.1]);
        q0 = quaternion(:, 1);
        q1 = quaternion(:, 2);
        q2 = quaternion(:, 3);
        q3 = quaternion(:, 4);
        val1 = 2 * ((q0 .* q1) + (q2 .* q3));
        val2 = q0.^2 - q1.^2 - q2.^2 + q3.^2;
        val3 = 2 .* ((q1 .* q2) + (q0 .* q3));
        val4 = 2 .* ((q1 .* q3) - (q0 .* q2));

        if (val1 < 1e-10) && (val2 < 1e-10)
            disp('Caution: Euler pitch angle close or at to 90 deg.')
        end
    end

    q0 = quaternion(:, 1);
    sq0 = q0.^2;
    q1 = quaternion(:, 2);
    sq1 = q1.^2;
    q2 = quaternion(:, 3);
    sq2 = q2.^2;
    q3 = quaternion(:, 4);
    sq3 = q3.^2;

    euler(:, 1) = atan2(2.0.*(q1 .* q2 + q3 .* q0), (sq1 - sq2 - sq3 + sq0));
    euler(:, 2) = asin(2.0.*(q2 .* q0 - q1 .* q3));
    euler(:, 3) = atan2(2.0.*(q2 .* q3 + q1 .* q0), (-sq1 - sq2 + sq3 + sq0));

    return
