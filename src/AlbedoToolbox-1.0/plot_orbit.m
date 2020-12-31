% PLOT_ORBIT Plots the orbit with REFL data background.
%
% plot_orbit(satpos,refl,type)
%
% where satpos is a matrix where each row is the position of the satellite
% in ECEF coordinates, refl is the REFL data to use as backgrpound, and type
% is either 'p' for flat plot or 's' for spherical plot.
%
% $Id: plot_orbit.m,v 1.5 2006/05/17 14:39:17 danji Exp $

function plot_orbit(satpos,refl,type)

if nargin < 3
  type = 'p';
end

if size(satpos,2) == 4
  warning('plot_orbit.m: satpos has dimension of 4. Sure time is not in there?');
end

plot_refl(refl,type,'no colorbar');
hold('on');
if strcmp(type,'p')
  [th,ph,r] = cart2sph(satpos(:,1),satpos(:,2),satpos(:,3));
  th = th.*180/pi;
  ph = ph.*180/pi;
  plot3(th,ph,ones(size(satpos,1),1)*100,'b.');
elseif strcmp(type,'s')
  for i=1:size(satpos,1)
    satpos(i,:) = satpos(i,:)./norm(satpos(i,:));
  end
  plot3(satpos(:,1),satpos(:,2),satpos(:,3),'b.');
else
  error('plot_orbit.m: Unknown plot type');
end

hold('off');
title('Orbit');

return
