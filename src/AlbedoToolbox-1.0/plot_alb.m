% PLOT_ALB Plot albedo data.
%
% h = plot_alp(albedo [,type])
%
% albedo is a matrix. Optional parameter type is 'p' for 3D
% plot or 's' for spherical. Default is 'p'.
%
% $Id: plot_alb.m,v 1.12 2006/05/17 14:39:17 danji Exp $

function h = plot_alb(albedo,type,cbar);

if nargin > 1 && strcmp(type,'s')
  if nargin > 2
    plot_alb_sphere(albedo,cbar);
  else
    plot_alb_sphere(albedo);
  end
	return;
end

[sy sx] = size(albedo);
dy = 180/sy;
dx = 360/sx;

lat = [-90+dy/2:dy:90-dy/2]';
lon = [-180+dx/2:dx:180-dx/2]';

hp = surf(lon,lat,albedo);
view(0,90);
axis([-180+dx 180-dx -90+dy 90-dy]);
shading('interp');
title('Albedo');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
zlabel('Irradiance [W/m^2]');
colormap(pink);

if ~(nargin > 2)
  % Add colorbar
  oldh = gca;
  cbhandle = colorbar;
  axes(cbhandle);
  ylabel('Irradiance [W/m^2]');
  axes(oldh);
end

% return handles, if requested
if nargout > 0
    h = gcf;
end

return
