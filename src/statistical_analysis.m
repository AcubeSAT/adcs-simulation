error_struct = load('instant_error_know');
error = error_struct.instant_error_know';
final_error = abs(error(:,:));

median_error(1) = median(final_error(1,:),'omitnan');
median_error(2) = median(final_error(2,:),'omitnan');
median_error(3) = median(final_error(3,:),'omitnan');

error_cov(1) = cov(final_error(1,:),'omitrows');
error_cov(2) = cov(final_error(2,:),'omitrows');
error_cov(3) = cov(final_error(3,:),'omitrows');

error_estimation(1) = median_error(1)+error_cov(1);
error_estimation(2) = median_error(2)+error_cov(2);
error_estimation(3) = median_error(3)+error_cov(3);

figure
subplot(3,1,1)
plot(final_error(1,:))
subplot(3,1,2)
plot(final_error(1,:))
subplot(3,1,3)
plot(final_error(3,:))
