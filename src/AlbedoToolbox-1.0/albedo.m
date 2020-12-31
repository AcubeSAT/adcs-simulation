% ALBEDO Calculation of albedo for a given satellite and Sun constellation
% and specified reflectivity data.
%
% a = albedo(sat, sun, refl [,type])
%
% sat and sun are the vectors to the Earth and Sun in Cartesian Earth Centered 
% Earth Fixed coordinates, respectively. refl is the reflectivity data to use 
% for albedo calculation. Type is 'p' for 3D plots and 's' for spherical. If
% unspecified no plots are generated.
%
% $Id: albedo.m,v 1.28 2006/05/17 14:39:16 danji Exp $

function a = albedo(sat,sun,refl,type);

CONST.EMR = 6371.01e3;
%CONST.AM0 = 1366.9;
CONST.AM0 = 1;
CONST.d2r = pi/180;

% Support row vectors

if length(sat) > size(sat,1)
  sat = sat';
end

if length(sun) > size(sun,1)
  sun = sun';
end

% Data size
[sy sx] = size(refl.data);

% Check for satellite attitude common error
if norm(sat) < CONST.EMR
  error('albedo.m: The satellite has crashed into Earth!');
end

% Spherical coordinates
[satsph(1) satsph(2) satsph(3)] = cart2sph(sat(1),sat(2),sat(3));
[sunsph(1) sunsph(2) sunsph(3)] = cart2sph(sun(1),sun(2),sun(3));

% Convert phi to polar angle
satsph(2) = pi/2 - satsph(2);
sunsph(2) = pi/2 - sunsph(2);

% REFL indices
[sun_i sun_j] = rad2idx(sunsph(1),sunsph(2),sy,sx);
[sat_i sat_j] = rad2idx(satsph(1),satsph(2),sy,sx);

% Visible elements
fov = earthfov(satsph,[sy sx]);
if nargin > 3
	figure(1);
	subplot(3,1,1);
	plot_refl(mask(refl.data,fov),type,'no colorbar');
	title('Satellite Field of View');
end

% Sunlit elements
vis = earthfov(sunsph,[sy sx]);
if nargin > 3
	subplot(3,1,2);
	plot_refl(mask(refl.data,vis),type,'no colorbar');
	title('Solar Field of View');
end

% Union
union = fov & vis;
if nargin > 3
	subplot(3,1,3);
	plot_refl(mask(refl.data,union),type,'no colorbar');
	title('Sunlit Satellite Field of View');
end

% Solar irradiance
irad = CONST.AM0; % [W/m^2]
index = 1;

% Loop through refl array and select sunlit and satellite visible cells
a = zeros(sy,sx);
grid = zeros(3,1);
for i=1:sy
    for j=1:sx
        if union(i,j)
					% Angle of incident solar irradiance
					phi_in = gridangle(i,j,sun_i,sun_j,sy,sx);
          % Account for numerical inaccuracies.
          if phi_in > pi/2
            phi_in = pi/2;
          end
					% Incident power					
					E_in = irad*cellarea(i,j,sy,sx)*cos(phi_in);
					% Distance to sat from grid
					[grid_theta grid_phi] = idx2rad(i,j,sy,sx);
					[grid(1) grid(2) grid(3)] = sph2cart(grid_theta,pi/2-grid_phi,CONST.EMR);
          %grid = [grid(1);grid(2);grid(3)];
					satdist = norm(sat-grid);
          % Angle to sat from grid
					phi_out = acos(((sat-grid)/satdist)'*grid/norm(grid));
          P_out = E_in*refl.data(i,j)*cos(phi_out)/(pi*satdist^2);
          % Store in array
				  a(i,j) = P_out;
          index = index + 1;
        end
    end
end

if nargin > 3
	figure(2);
	plot_alb(a,type);
end

return
