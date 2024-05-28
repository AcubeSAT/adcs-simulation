%% Sun field of view plots
%
% Inputs:
%  quaternion_ECI_to_Body   - The quaternion from ECI to Body frame
%  Time                     - the time that we want to find the field of view
%  sun_pos_eci              - the position of the sun in the ECI Frame
%
% Both quaternion and the sun vector need to be matrices, derived from a
% whole orbit (or more), to produce the plot for this interval.
%
%
% The purpose of this function is to calculate the angle between a certain
% plain of the satellite and the sun. We can select the plain that we want
% to examine by choosing the proper euler angles in the function
% angles_to_quaternion. Then, by turning the quaternion from body to the
% plane that we want, we find the new frame between the sun and the
% satellite's plane. The sun vector is calculated in this new frame with
% the use of quaternions. This means that we calculate a new quaternion,
% sun_pos_sensor, where we keep only the final 3 elements, which are the
% three components of the vector, x,y,z.
% The angle between this plain and the sun is calculated as follows:
% cos(a) = u_x/norm(u) where a is the angle between x axis and u
% In our case, a is equal to the angle between the satellite's plane and
% the sun and the vector u is the sun vector in this new frame. Note that
% we only want the projection of the vector in the axis perpendicular to
% the plane, hence we only keep the first component of the vector, or the
% second component of the quaternion sun_pos_sensor. Finally we transform
% this angle from rads to degrees.


function sun_fov_plots(quaternion_ECI_to_Body, Time, sun_pos_eci)

    quat_sensor = zeros(4, length(Time)+1);
    sun_pos_sensor = zeros(4, length(Time)+1);

    q_rotation = angles_to_quaternion(pi/2, 0, 0);

    for i = 1:length(quaternion_ECI_to_Body(1, :))
        quat_sensor(:, i) = quatProd(quaternion_ECI_to_Body(1:4, i), q_rotation);
    end


    for i = 1:length(Time)
        sun_pos_sensor(:, i) = quatProd(quatconj(quat_sensor(:, i)'), quatProd([0; sun_pos_eci(:, i)], quat_sensor(:, i)));

        sun_pos_sensor(:, i) = sun_pos_sensor(:, i) / norm(sun_pos_sensor(:, i));
    end


    azimuth = zeros(1, length(Time));
    elevation = zeros(1, length(Time));
    fov_angle = zeros(1, length(Time));

    for i = 1:length(Time)
        [azimuth(i), elevation(i)] = cart2sph(sun_pos_sensor(2, i), sun_pos_sensor(3, i), sun_pos_sensor(4, i));
        azimuth = azimuth * 180 / pi;
        elevation = elevation * 180 / pi;
        cosine_fov_angle = sun_pos_sensor(2, i);
        fov_angle(i) = acos(cosine_fov_angle);
    end
    fov_angle = fov_angle * 180 / pi;

    figure();
    plot(Time, fov_angle, 'LineWidth', 2.0, 'Color', 'blue')
    title('Angle between sun vector & sun sensor');
    xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 15);
    ylabel('[$angle$]', 'interpreter', 'latex', 'fontsize', 17);

end