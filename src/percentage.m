final = length(instant_error_perform);
indices = [1:15962 (29525+5000):72597 (86242+5000):129072 (142960+5000):186050 (199677+5000):242583 (256394+5000):299561 (313112+5000):356170 (369829+5000):412339 (426546+5000):469430 (483263+5000):526519 (539980+5000):final]; %--> 6pm
% indices = [1:22092 (43496):78813 (100217):135533 (156937):192189 (213657):248589 (270377):305525 (327097):332701]; --> 11pm
% indices = [1:22611 (42025):80668 (100081):138725 (158136):196781 (216192):254838 (274248):312894 332304:332701]; --> 600km
% indices = [1:20721 (40676):76918 (96873):133115 (153070):189313 (209267):245510 (265463):301707 (321660):332701]; --> 450km
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
y = abs(angles(:,2)) <= 40;

total1 = 0;
total2 = 0;
total3 = 0;

for i=1:length(y)
    total1 = total1 + x(i);
    total2 = total2 + y(i);
    total3 = total3 + y(i)*x(i);
    
end

per1 = total1/length(y) %X < 20
per2 = total2/length(y) %Y < 20
per3 = total3/length(y) % Z < 20

error_estimation = zeros(2,2);

alpha = 0.05;

error_estimation(1,:) = bootci(100,{@median,angles(:,1)},'alpha',alpha,'type','percentile');
error_estimation(2,:) = bootci(100,{@median,angles(:,2)},'alpha',alpha,'type','percentile')

