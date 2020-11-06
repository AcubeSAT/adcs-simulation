function sun_fov_plots(x_data,Time,sun_pos_eci)

quat_sensor=zeros(4,length(Time)+1);
sun_pos_sensor = zeros(4,length(Time)+1);


% Calculate sun sensor reference frame    Reflection along YZ plane = -w -x
quat_sensor(1,:) = -x_data(1,:);
quat_sensor(2,:) = -x_data(2,:);
quat_sensor(3,:) = x_data(3,:);
quat_sensor(4,:) = x_data(4,:);

for i=1:length(Time)   % Calculate Sun in sun sensor reference frame
    sun_pos_sensor(:,i) = quatProd(quatconj(quat_sensor(:,i)'),quatProd( [0;sun_pos_eci(:,i)] , quat_sensor(:,i)));

    sun_pos_sensor(:,i) =  sun_pos_sensor(:,i)/norm(sun_pos_sensor(:,i)); % Normalize
end


quat_fov=zeros(4,length(Time));
fov_angle=zeros(1,length(Time));
for i=1:length(Time)
quat_fov(:,i)=quatProd(quatconj(quat_sensor(:,i)'),sun_pos_sensor(:,i)); % Calculate FOV using Sun_body
fov_angle(1,i) = 2*asin(norm(quat_fov(2:4,i)))*180/pi; 

end

figure();
plot(Time,fov_angle(1,:),'LineWidth',2.0, 'Color','blue')
title('Angle between sun vector & sun sensor');
xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',15);
ylabel(['$angle$'], 'interpreter','latex', 'fontsize',17);
end