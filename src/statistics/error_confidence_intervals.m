%% ======================================================================= %%
% This script calculates the confidence intervals of median absolute error
% on each axis, before and after filtering out the eclipse periods.
% Also, it plots the absolute performace error.
% It is to be executed only after Nadir Pointing simulation.
%
%   Input  : Performance/knowledge error matrix
%   Output : 95% confidence intervals of median absolute error on each axis
% ======================================================================== %

%% Exclude periods of eclipse
%offset = 5000;
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

error = instant_error_perform';
filtered_error = instant_error_perform(indices, :)';

%% Main script
% error = instant_error_perform';  % Input matrix
final_error = abs(error);
final_error = final_error(:, all(~isnan(final_error)));

final_filtered_error = abs(filtered_error);
final_filtered_error = final_filtered_error(:, all(~isnan(final_filtered_error)));

error_estimation = zeros(3, 2);
filtered_error_estimation = zeros(3, 2);
alpha = 0.05;

error_estimation(1, :) = bootci(100, {@median, final_error(1, :)}, 'alpha', alpha, 'type', 'percentile')';
error_estimation(2, :) = bootci(100, {@median, final_error(2, :)}, 'alpha', alpha, 'type', 'percentile')';
error_estimation(3, :) = bootci(100, {@median, final_error(3, :)}, 'alpha', alpha, 'type', 'percentile')';

error_estimation

filtered_error_estimation(1, :) = bootci(100, {@median, final_filtered_error(1, :)}, 'alpha', alpha, 'type', 'percentile')';
filtered_error_estimation(2, :) = bootci(100, {@median, final_filtered_error(2, :)}, 'alpha', alpha, 'type', 'percentile')';
filtered_error_estimation(3, :) = bootci(100, {@median, final_filtered_error(3, :)}, 'alpha', alpha, 'type', 'percentile')';

filtered_error_estimation

figure
sgtitle('Absolute performance error on each axis', 'interpreter', 'latex', 'fontsize', 18)
subplot(3, 1, 1)
plot(final_error(1, :))
hold on;
plot(50*eclipse, 'LineWidth', 1.5, 'Color', 'green')
hold on
yline(20, 'LineWidth', 1.5, 'Color', 'red')
ylabel('Z-axis', 'interpreter', 'latex', 'fontsize', 14)
subplot(3, 1, 2)
plot(final_error(2, :))
hold on;
plot(25*eclipse, 'LineWidth', 1.5, 'Color', 'green')
hold on
yline(20, 'LineWidth', 1.5, 'Color', 'red')
ylabel('Y-axis', 'interpreter', 'latex', 'fontsize', 14)
subplot(3, 1, 3)
plot(final_error(3, :))
hold on;
plot(50*eclipse, 'LineWidth', 1.5, 'Color', 'green')
hold on
yline(20, 'LineWidth', 1.5, 'Color', 'red')
ylabel('X-axis', 'interpreter', 'latex', 'fontsize', 14)
