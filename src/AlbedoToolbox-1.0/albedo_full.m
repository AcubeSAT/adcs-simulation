% ALBEDO_FULL Calculate albedo array for 100% sunlight (zenith) at all
% satellite positions at a given altitude. The resulting albedo is the
% maximum albedo for all satellite positions.
%
% a = albedo_full(altitude, refl)
%
% altitude is the altitude of the satellite in meters. refl is the
% reflectivity data.
%
% $Id: albedo_full.m,v 1.7 2006/05/17 14:39:16 danji Exp $

function a = albedo_full(altitude,refl);

CONST.EMR = 6371.01e3;
CONST.AU = 149598000000;

% Data size
[sy sx] = size(refl.data);
sat = zeros(3,1);
sun = zeros(3,1);

h = waitbar(0,'Calculating full albedo matrix...');
for i = 1:sy
	for j = 1:sx
		[sat_theta sat_phi] = idx2rad(i,j,sy,sx);
		[sat(1) sat(2) sat(3)] = sph2cart(sat_theta,sat_phi-pi/2,CONST.EMR+altitude);
		[sun(1) sun(2) sun(3)] = sph2cart(sat_theta,sat_phi-pi/2,CONST.AU);
		a(i,j) = sum(sum(albedo(sat,sun,refl)));
		waitbar(i*j/(sy*sx),h);
	end
end

close(h);

return
