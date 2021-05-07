%Input: Performance/knowledge error matrix
%Output: 95% confidence intervals of median absolute error on each axis

%% The following is for excluding periods of eclipse
%error_struct = load('performance_error1.mat');
% indices = [1:15962 (29525):72540 (86242):129072 (142960):166351 186050 (199677):242583 (256394):299561 (313112):332701]; --> 6pm
% indices = [1:22092 (43496):78813 (100217):135533 (156937):192189 (213657):248589 (270377):305525 (327097):332701]; --> 11pm
% indices = [1:22611 (42025):80668 (100081):138725 (158136):196781 (216192):254838 (274248):312894 332304:332701]; --> 600km
% indices = [1:20721 (40676):76918 (96873):133115 (153070):189313 (209267):245510 (265463):301707 (321660):332701]; --> 450km
% indices = [1:23415 (44639):81471 (102695):139527 (160751):197584 (218808):255640 (276864):length(instant_error_perform)]; % --> 600km 11pm

error = instant_error_perform(indices,:)';

%% Main script
% error = instant_error_perform';  % Input matrix
final_error = abs(error);
final_error = final_error(:,all(~isnan(final_error)));
error_estimation = zeros(3,2);
alpha = 0.05;

error_estimation(1,:) = bootci(100,{@median,final_error(1,:)},'alpha',alpha,'type','percentile')';
error_estimation(2,:) = bootci(100,{@median,final_error(2,:)},'alpha',alpha,'type','percentile')';
error_estimation(3,:) = bootci(100,{@median,final_error(3,:)},'alpha',alpha,'type','percentile')';

error_estimation

figure
subplot(3,1,1)
plot(final_error(1,:))
subplot(3,1,2)
plot(final_error(1,:))
subplot(3,1,3)
plot(final_error(3,:))
