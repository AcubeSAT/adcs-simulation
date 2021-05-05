
anles = zeros(length(instant_error_perform),2);

for i = 1:length(instant_error_perform)
q=eul2quat(deg2rad(instant_error_perform(i,:)));
R_OB = quat2dcm(q);
R_BO =R_OB';
keraia_vector = (R_BO*[1 0 0]');
[azimuth,elevation,r] = cart2sph(keraia_vector(1),keraia_vector(2),keraia_vector(3));

angles(i,1) = abs(rad2deg(azimuth));
angles(i,2) = abs(rad2deg(elevation));
end
figure()
plot(angles(:,1))
figure()
plot(angles(:,2))

