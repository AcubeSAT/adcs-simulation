%% CSS Noise Function
%
% Inputs:
%  sun_eci      - Sun in the ECI frame
%  q_eci_body   - the quaternion from ECI to Body frame
%  xsat_eci     - the position of the Satellite in the ECI frame
%  albedo_perc  - the percentage of the sunlight diffused from earth that
%                 interacts with the satellite
%  lambda       - the poisson parameter used to add noise
%
% Outputs:
%  total_sun_vector - the total sun vector, i.e. including both sun, albedo,
%                     and error measurement in the body frame
%
%
% The function finds where each of the 6 Coarse sun sensors look at the
% giving time by turning the Body frame to an angle according to the
% position of the sensor, and then calculates how much albedo each sensor
% receives and adds this as a noise in the measurement. All albedo radiation
% is asumed to come from Nadir. Furthermore, in all 6 measurements,
% a Poisson noise is added, since CSS are not 100% accurate.
% The measurements of the CSS are 6 currents, from which we can calculate
% the total sun vector in the Body Frame. If one of the currents is negative,
% it means that this sensor does not see the sun.
% For more details, you can refer to DDJF_AOCS file.

function total_sun_vector = css_noise(sun_eci, q_eci_body, xsat_eci, albedo_perc, lambda)

    sun_eci = sun_eci / norm(sun_eci);
    sun_body = rotate_vector(q_eci_body, sun_eci);

    xsat_eci = xsat_eci / norm(xsat_eci);
    xsat_body = rotate_vector(q_eci_body, xsat_eci);
    nadir = -xsat_body;

    %frame_of_css_1 = roty(0); the first frame is the same with the Body!
    % frame_of_css_2 = roty(90);
    % frame_of_css_3 = roty(-90);
    % frame_of_css_4 = roty(180);
    % frame_of_css_5 = rotz(90);
    % frame_of_css_6 = rotz(-90);

    % Use custom functions instead of toolbox
    frame_of_css_2 = rot_y(90);
    frame_of_css_3 = rot_y(-90);
    frame_of_css_4 = rot_y(180);
    frame_of_css_5 = rot_z(90);
    frame_of_css_6 = rot_z(-90);

    %q1 = dcm2quat(f1);
    q2 = dcm2quat(frame_of_css_2);
    q3 = dcm2quat(frame_of_css_3);
    q4 = dcm2quat(frame_of_css_4);
    q5 = dcm2quat(frame_of_css_5);
    q6 = dcm2quat(frame_of_css_6);

    final_sun(:, 1) = sun_body;
    final_sun(:, 2) = rotate_vector(q2', sun_body);
    final_sun(:, 3) = rotate_vector(q3', sun_body);
    final_sun(:, 4) = rotate_vector(q4', sun_body);
    final_sun(:, 5) = rotate_vector(q5', sun_body);
    final_sun(:, 6) = rotate_vector(q6', sun_body);

    final_nadir(:, 1) = nadir;
    final_nadir(:, 2) = rotate_vector(q2', nadir);
    final_nadir(:, 3) = rotate_vector(q3', nadir);
    final_nadir(:, 4) = rotate_vector(q4', nadir);
    final_nadir(:, 5) = rotate_vector(q5', nadir);
    final_nadir(:, 6) = rotate_vector(q6', nadir);


    current_sun = zeros(6, 1);
    current_albedo = zeros(6, 1);
    total_current = zeros(6, 1);

    for i = 1:6
        current_sun(i) = dot(final_sun(:, i), [1; 0; 0]);
        sign = randi([0, 1]);
        if sign == 0
            sign = -1;
        end
        if current_sun(i) < 0
            current_sun(i) = 0;
        end
        current_albedo(i) = albedo_perc * dot(final_nadir(:, i), [1; 0; 0]);
        if current_albedo(i) < 0
            current_albedo(i) = 0;
        end
        total_current(i) = current_sun(i) + current_albedo(i);
        total_current(i) = total_current(i) + total_current(i) * sign * 0.01 * poissrnd(lambda);
    end

    total_sun_vector = [total_current(1) - total_current(4), total_current(6) - total_current(5), total_current(2) - total_current(3)];
    total_sun_vector = (total_sun_vector / norm(total_sun_vector))'; % Normalized sun vector in body frame

end
