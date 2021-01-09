% DOUBLE_PLOT_ALB Plot albedo data on subplot with dual view.
%
% double_plot_alb(a)
%
% where a is an albedo matrixx.
%
% $Id: double_plot_alb.m,v 1.3 2006/02/23 08:31:33 danji Exp $

function double_plot_alb(a);

subplot(2,1,1);
plot_alb(a);

subplot(2,1,2);
[sy sx] = size(a);
dy = 180/sy;
dx = 360/sx;

lat = [-90+dy/2:dy:90-dy/2]';
lon = [-180+dx/2:dx:180-dx/2]';

hp = surf(lon,lat,a);
axis([-180+dx 180-dx -90+dy 90-dy]);
shading('interp');
title('Albedo');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
zlabel('Irradiance [W/m^2]');
view(-34,32);
colormap(pink);

return
