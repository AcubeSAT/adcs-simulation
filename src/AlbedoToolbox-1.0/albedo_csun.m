% ALBEDO_CSUN Calculate albedo array for constant sunlight (instantanious
% albedo) for all satellite positions at given altitude.
%
% a = albedo_csun(altitude, sunsph, refl [,type])
%
% altitude is the altitude of the albedo calculations in meters. sunsph is the
% position of the Sun in ECEF spherical coordinates. refl is the
% reflectivity data. Type is 'p' for 3D plot and 's' for spherical. If
% unspecified no plot is generated.
%
% $Id: albedo_csun.m,v 1.9 2006/05/17 14:39:16 danji Exp $

function a = albedo_csun(altitude,sunsph,refl,type);

CONST.EMR = 6371.01e3;

% Data size
[sy sx] = size(refl.data);

if nargin > 3
  % Sunlit elements
  vis = earthfov(sunsph,refl,type);
  plot_refl(mask(refl.data,vis),type);
end

sat = zeros(3,1);
sun = zeros(3,1);

[sun(1) sun(2) sun(3)] = sph2cart(sunsph(1),sunsph(2)-pi/2,sunsph(3));

h = waitbar(0,'Calculating instantanious constant sun albedo matrix...');

for i = 1:sy
	for j = 1:sx
		[sat_theta sat_phi] = idx2rad(i,j,sy,sx);
		[sat(1) sat(2) sat(3)] = sph2cart(sat_theta,sat_phi-pi/2,CONST.EMR+altitude);
		a(i,j) = sum(sum(albedo(sat,sun,refl)));
		waitbar(i*j/(sy*sx),h);
  end
end

close(h);

return
