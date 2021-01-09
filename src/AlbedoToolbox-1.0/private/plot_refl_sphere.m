% PLOT_REFL_SPHERE Plot REFL data on a sphere.
%
% plot_refl_sphere(refl_parm)
%
% refl_parm is a struct containing a 180x288 matrix in a
% field named data or a 180x288 matrix. Add second argument to omit
% colorbar.
%
% $Id: plot_refl_sphere.m,v 1.19 2006/05/17 14:39:18 danji Exp $

function h = plot_refl_sphere(refl_parm,cbar)

if ~isfield(refl_parm,'data')
  refl.data = refl_parm;
else
  refl = refl_parm;
end

d2r = pi/180;

[sy sx] = size(refl.data);
dy = 180/sy;
dx = 360/sx;

phi = [-90+dy/2:dy:90-dy/2]'.*d2r;
theta = [-180+dx/2:dx:180-dx/2]'.*d2r;

X = cos(theta)*cos(phi)';
Y = sin(theta)*cos(phi)';
Z = ones(sx+1,1)*sin(phi)';

X = [X;X(1,:)];
Y = [Y;Y(1,:)];

hp = surf(X,Y,Z,zeros(sx+1,sy));
hold on;
hp = surf(X,Y,Z,[refl.data'.*100;refl.data(:,1)'.*100]);
hold off;
shading('interp');
title('TOMS Reflectivity Data');
view(127.5,16);
colormap(pink);
% Remove useless ticks
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])

if ~(nargin > 1)
  % Add colorbar
  oldh = gca;
  cbhandle = colorbar;
  axes(cbhandle);
  ylabel('Reflectivity [%]');
  axes(oldh);
end

if nargout > 0
    % Return handle if requested
    h = gcf;
end

return
