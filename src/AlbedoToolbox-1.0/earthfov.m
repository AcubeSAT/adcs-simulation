% EARTHFOV Field of view on earth by spherical coordinates.
%
% result = earthfov(satsph, refl [, type])
%
% satsph is the vector to the satellite in ECEF frame and spherical
% coordinates. refl is the reflectivity data used for plotting the fov.
% type is 'p' for 3D plot and 's' for spherical. If refl or type are
% unspecified no plot is generated. If refl is not specified, refl must be 
% a vector of the latitudal (refl(1)) and longitudal (refl(2)) resolution 
% of the output data.
%
% $Id: earthfov.m,v 1.15 2006/02/23 08:31:33 danji Exp $

function result = earthfov(satsph,refl,type);

CONST.EMR = 6371.01e3;
CONST.d2r = pi/180;

% LEO shortcut
if satsph(3) < CONST.EMR
	satsph(3) = satsph(3) + CONST.EMR;
end

% Circle value
OUTVALUE = 1;

% Discretization parameters
if isfield(refl,'data')
  [sy sx] = size(refl.data);
else 
  sy = refl(1);
  sx = refl(2);
end
dx = 2*pi/sx;
dy = pi/sy;
result = zeros(sy,sx);

% Small Circle Center
theta0 = satsph(1);
phi0 = satsph(2);

% FOV on earth
rho = acos(CONST.EMR/satsph(3));

for i = 1:sy
	for j = 1:sx        
		[theta, phi] = idx2rad(i,j,sy,sx);
		% Radial distance
		rd = acos(sin(phi0)*sin(phi)*cos(theta0-theta)+cos(phi0)*cos(phi));
		if rd <= rho
			result(i,j) = OUTVALUE;
		end
	end
end

if nargin > 1 && isfield(refl,'data')
  if nargin > 2
    plot_refl(mask(refl.data,result),type,'no colorbar');
  else
    plot_refl(mask(refl.data,result),'p','no colorbar');
  end
  title('Field of View');
end

return
