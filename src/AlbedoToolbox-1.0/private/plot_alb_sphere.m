% PLOT_ALB_SPHERE Plot albedo data on a sphere.
%
% plot_alb_sphere(albedo)
%
% albedo is a 180x288 matrix.
%
% $Id: plot_alb_sphere.m,v 1.10 2006/05/17 14:39:18 danji Exp $

function h = plot_alb_sphere(albedo)

d2r = pi/180;

[sy sx] = size(albedo);
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
hp = surf(X,Y,Z,[albedo';albedo(:,1)']);
hold off;
shading('interp');
title('Albedo');
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
  ylabel('Irradiance [W/m^2]');
  axes(oldh);
end

if nargout > 0
    % Return handle if requested
    h = gcf;
end

return
