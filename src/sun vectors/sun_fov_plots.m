function sun_fov_plots(x_data,Time,sun_pos_eci)

quat_sensor=zeros(4,length(Time)+1);
sun_pos_sensor = zeros(4,length(Time)+1);

q_rotation = angles_to_quaternion(pi/2,0,0);
for i=1:length(x_data(1,:))
    quat_sensor(:,i) = quatProd(x_data(1:4,i), q_rotation);
end


for i=1:length(Time)   % Calculate Sun in sun sensor reference frame
    sun_pos_sensor(:,i) = quatProd(quatconj(quat_sensor(:,i)'),quatProd( [0;sun_pos_eci(:,i)] , quat_sensor(:,i)));

    sun_pos_sensor(:,i) =  sun_pos_sensor(:,i)/norm(sun_pos_sensor(:,i)); % Normalize
end



azimuth = zeros(1,length(Time));
elevation = zeros(1,length(Time));
fov_angle = zeros(1,length(Time));
for i=1:length(Time)
    [azimuth(i),elevation(i)] = cart2sph(sun_pos_sensor(2,i),sun_pos_sensor(3,i),sun_pos_sensor(4,i));
    azimuth = azimuth*180/pi;
    elevation = elevation*180/pi;
    cosine_fov_angle = sun_pos_sensor(2,i); %cos(a) = u_x/norm(u) where a is the angle between x axis and u
    fov_angle(i) = acos(cosine_fov_angle); 
end
fov_angle = fov_angle*180/pi;

figure();
plot(Time,fov_angle,'LineWidth',2.0, 'Color','blue')
title('Angle between sun vector & sun sensor');
xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',15);
ylabel(['$angle$'], 'interpreter','latex', 'fontsize',17);
end