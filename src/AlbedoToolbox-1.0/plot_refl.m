% PLOT_REFL Plot REFL data.
%
% plot_refl(refl_parm [, type [, colorbar]])
%
% refl_parm is a struct containing a matrix in a field named data or a
% matrix itself. Optional parameter type is 'p' for 3D plot or 's' for
% spherical. Default is 'p'. To leave out colorbar, add third argument
% to function call.
%
% $Id: plot_refl.m,v 1.22 2006/05/17 14:39:17 danji Exp $

function h = plot_refl(refl_parm,type,cbar);

if ~isfield(refl_parm,'data')
  refl.data = refl_parm;
else
  refl = refl_parm;
end
  
if nargin > 1 && strcmp(type,'s')
  if nargin > 2
    plot_refl_sphere(refl,cbar);
  else
    plot_refl_sphere(refl);
  end
	return;
end

[sy sx] = size(refl.data);
dy = 180/sy;
dx = 360/sx;

lat = [-90+dy/2:dy:90-dy/2]';
lon = [-180+dx/2:dx:180-dx/2]';

hp = surf(lon,lat,refl.data.*100);
view(0,90);
axis([-180+dx 180-dx -90+dy 90-dy]);
shading('interp');
title('TOMS Reflectivity Data');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
zlabel('Reflectivity [%]');
colormap(pink);

if ~(nargin > 2)
  % Add colorbar
  oldh = gca;
  cbhandle = colorbar;
  axes(cbhandle);
  ylabel('Reflectivity [%]');
  axes(oldh);
end

% return handles, if requested
if nargout > 0
    h = gcf;
end

return
