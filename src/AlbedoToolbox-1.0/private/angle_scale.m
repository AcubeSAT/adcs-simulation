% ANGLE_SCALE See index to angle conversion.
%
% angle_scale(sy,sx)
%
% where sy and sx are the dimensions of the latitude and longitude scales.
%
% $Id: angle_scale.m,v 1.2 2006/05/17 14:39:17 danji Exp $

function angle_scale(sy,sx)

r2d = 180/pi;

n=1000;
idx = zeros(n+1,2);

for i=1:n+1
  [idx(i,1) idx(i,2)] = rad2idx((i-1)*2*pi/n-pi,(i-1)*pi/n,sy,sx);
end

subplot(2,1,1);
plot([0:pi/n:pi]*r2d,idx(:,1));
title('Elevation map');
subplot(2,1,2);
plot([-pi:2*pi/n:pi]*r2d,idx(:,2));
title('Azimuth map');

return
