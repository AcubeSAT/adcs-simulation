final = length(instant_error_perform);
offset = 0;
timestamps= [];
e =1;
for i = 1:length(eclipse)
    if eclipse(i) == e
        timestamps = [timestamps i];
        if e == 1
            e = 0;
        else
            e = 1;
        end
    end
end
if mod(length(timestamps),2) == 0
    timestamps = [timestamps 0];
end
indices = 1:timestamps(1);


for i = 2:2:length(timestamps)
    if timestamps(i+1) ~=0
        indices = [indices (timestamps(i)+offset):timestamps(i+1)];
    else
        indices = [indices (timestamps(i)+offset):length(instant_error_perform)];
    end
end

error = instant_error_perform(indices,:);
angles = zeros(length(error),2);

for i = 1:length(error)
    
q=eul2quat(deg2rad(error(i,:)));
R_OB = quat2dcm(q);
R_BO =R_OB';
antenna_vector = (R_BO*[1 0 0]');

[azimuth,elevation,~] = cart2sph(antenna_vector(1), antenna_vector(2), antenna_vector(3));

angles(i,1) = abs(azimuth*180/pi);
angles(i,2) = abs(elevation*180/pi);

end

x = abs(angles(:,1)) <= 20;
y = abs(angles(:,2)) <= 20;
z = zeros(1,length(angles(:,1)));
total1 = 0;
total2 = 0;
total3 = 0;

for i=1:length(y)
    total1 = total1 + x(i);
    total2 = total2 + y(i);
    total3 = total3 + y(i)*x(i);
    if y(i)*x(i) == 1
        z(i) = 1;
    end
end

per1 = total1/length(y) %X < 20
per2 = total2/length(y) %Y < 20
per3 = total3/length(y) % Z < 20

error_estimation = zeros(2,2);

z_last = z(1);
i_last = 0;

for i = 1 :length(z)
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




% figure()
% plot(z,'LineWidth',1.5, 'Color','blue')
% hold on;
% plot(eclipse, 'LineWidth',1.5, 'Color','magenta')
% ylim([0 3])


alpha = 0.05;

error_estimation(1,:) = bootci(100,{@median,angles(:,1)},'alpha',alpha,'type','percentile');
error_estimation(2,:) = bootci(100,{@median,angles(:,2)},'alpha',alpha,'type','percentile')


total = 0;
for i = 1:length(z)
if z(i)==1
total = total +1;
end
end
total/length(z)











