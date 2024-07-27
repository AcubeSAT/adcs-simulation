%% ======================================================================= %%
% This script calculates the percentage of timesteps where the azimuth
% and elevation angle are within range, as well as the confidence intervals,
% before and after filtering out the eclipse periods.
% Also, it plots
% It is to be executed only after Nadir Pointing simulation.
%
%   Input  : Performance/knowledge error matrix
%   Output : Percentage of timesteps where azimuth and elevation angle are
%            less than 20 degrees
%            95% confidence intervals of azimuth and elevation angle
% ======================================================================== %

%% Exclude periods of eclipse
%final = length(instant_error_perform);
offset = 0;
timestamps = [];
e = 1;
for i = 1:length(eclipse)
    if eclipse(i) == e
        timestamps = [timestamps, i];
        if e == 1
            e = 0;
        else
            e = 1;
        end
    end
end
if mod(length(timestamps), 2) == 0
    timestamps = [timestamps, 0];
end
indices = 1:timestamps(1);


for i = 2:2:length(timestamps)
    if timestamps(i+1) ~= 0
        indices = [indices, (timestamps(i) + offset):timestamps(i+1)];
    else
        indices = [indices, (timestamps(i) + offset):length(instant_error_perform)];
    end
end

error = instant_error_perform;
filtered_error = instant_error_perform(indices, :);

%% Calculate azimuth and elevation angles for patch antenna pointing
angles = zeros(length(error), 2);
filtered_angles = zeros(length(filtered_error), 2);

for i = 1:length(error)
    q = eul2quat(deg2rad(error(i, :)));
    R_OB = quat2dcm(q);
    R_BO = R_OB';
    antenna_vector = (R_BO * [1, 0, 0]');

    [azimuth, elevation, ~] = cart2sph(antenna_vector(1), antenna_vector(2), antenna_vector(3));

    angles(i, 1) = abs(azimuth*180/pi);
    angles(i, 2) = abs(elevation*180/pi);
end

for i = 1:length(filtered_error)
    q = eul2quat(deg2rad(filtered_error(i, :)));
    R_OB = quat2dcm(q);
    R_BO = R_OB';
    antenna_vector = (R_BO * [1, 0, 0]');

    [azimuth, elevation, ~] = cart2sph(antenna_vector(1), antenna_vector(2), antenna_vector(3));

    filtered_angles(i, 1) = abs(azimuth*180/pi);
    filtered_angles(i, 2) = abs(elevation*180/pi);
end

%% Main script for raw data
% Calculate if azimuth and elevation angles are less than or equal to 20 degrees
x = abs(angles(:, 1)) <= 20; % azimuth
y = abs(angles(:, 2)) <= 20; % elevation
z = zeros(1, length(angles(:, 1))); % both azimuth & elevation

total1 = 0;
total2 = 0;
total3 = 0;

% Calculate percentage of occurences where azimuth and elevation angles are within range
for i = 1:length(y)
    total1 = total1 + x(i);
    total2 = total2 + y(i);
    total3 = total3 + y(i) * x(i);
    if y(i) * x(i) == 1
        z(i) = 1;
    end
end

per1 = total1 / length(y) %X < 20
per2 = total2 / length(y) %Y < 20
per3 = total3 / length(y) %Z < 20

error_estimation = zeros(2, 2);

% Filter out short intervals of change in the array z
z_last = z(1);
i_last = 0;

for i = 1:length(z)
    if z(i) ~= z_last
        dif = i - i_last;
        z_last = z(i);
        if dif < 900 && z_last == 0
            for j = i:-1:i_last
                z(j) = 0;
            end
        end
        i_last = i;
    end
end


figure()
sgtitle('Azimuth and Elevation angle threshold', 'interpreter', 'latex', 'fontsize', 18)
plot(z, 'LineWidth', 1.5, 'Color', 'blue')
hold on;
plot(eclipse, 'LineWidth', 1.5, 'Color', 'magenta')
ylim([0, 3])


alpha = 0.05;

error_estimation(1, :) = bootci(100, {@median, angles(:, 1)}, 'alpha', alpha, 'type', 'percentile');
error_estimation(2, :) = bootci(100, {@median, angles(:, 2)}, 'alpha', alpha, 'type', 'percentile');
error_estimation

total = 0;
for i = 1:length(z)
    if z(i) == 1
        total = total + 1;
    end
end
total / length(z)

%% Main script for filtered data
% Calculate if azimuth and elevation angles are less than or equal to 20 degrees
filtered_x = abs(filtered_angles(:, 1)) <= 20; % azimuth
filtered_y = abs(filtered_angles(:, 2)) <= 20; % elevation
filtered_z = zeros(1, length(filtered_angles(:, 1))); % both azimuth & elevation

filtered_total1 = 0;
filtered_total2 = 0;
filtered_total3 = 0;

% Calculate percentage of occurences where azimuth and elevation angles are within range
for i = 1:length(filtered_y)
    filtered_total1 = filtered_total1 + filtered_x(i);
    filtered_total2 = filtered_total2 + filtered_y(i);
    filtered_total3 = filtered_total3 + filtered_y(i) * filtered_x(i);
    if filtered_y(i) * filtered_x(i) == 1
        filtered_z(i) = 1;
    end
end

filtered_per1 = filtered_total1 / length(filtered_y) %X < 20
filtered_per2 = filtered_total2 / length(filtered_y) %Y < 20
filtered_per3 = filtered_total3 / length(filtered_y) %Z < 20

filtered_error_estimation = zeros(2, 2);

% Filter out short intervals of change in the array z
z_last = filtered_z(1);
i_last = 0;

for i = 1:length(filtered_z)
    if z(i) ~= z_last
        dif = i - i_last;
        z_last = z(i);
        if dif < 900 && z_last == 0
            for j = i:-1:i_last
                z(j) = 0;
            end
        end
        i_last = i;
    end
end

alpha = 0.05;

filtered_error_estimation(1, :) = bootci(100, {@median, filtered_angles(:, 1)}, 'alpha', alpha, 'type', 'percentile');
filtered_error_estimation(2, :) = bootci(100, {@median, filtered_angles(:, 2)}, 'alpha', alpha, 'type', 'percentile');
filtered_error_estimation

filtered_total = 0;
for i = 1:length(filtered_z)
    if z(i) == 1
        filtered_total = filtered_total + 1;
    end
end
filtered_total / length(filtered_z)
