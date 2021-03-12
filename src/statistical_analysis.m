%error_struct = load('performance_error1.mat');
%indices = [1:15962 (29525+5000):72540 (86242+5000):129072 (142960+5000):166351 (199677):242583 (256394):299561 (313112):332701];
%error = instant_error_perform(indices,:)';
error = instant_error_perform';
final_error = abs(error(:,:));
final_error = final_error(:,all(~isnan(final_error)));
error_estimation = zeros(3,2);
alpha = 0.05;

error_estimation(1,:) = bootci(100,{@median,final_error(1,:)},'alpha',alpha,'type','percentile')';
error_estimation(2,:) = bootci(100,{@median,final_error(2,:)},'alpha',alpha,'type','percentile')';
error_estimation(3,:) = bootci(100,{@median,final_error(3,:)},'alpha',alpha,'type','percentile')';
figure
subplot(3,1,1)
plot(final_error(1,:))
subplot(3,1,2)
plot(final_error(1,:))
subplot(3,1,3)
plot(final_error(3,:))
