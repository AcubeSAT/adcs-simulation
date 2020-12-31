% ALBEDO_ALTITUDE Calculate subsolar albedo between two altitudes (above
% earth) at specified position.
%
% result = albedo_altitude(az, pa, alt1, alt2, refl, n [,type])
%
% az and pa is the azimuth and polar angle of the position in ECEF frame.
% Albedo is calculated from alt1 to alt2 with n intermediate calculations
% (defaults to 10). Note that altitudes must be given in meters. Type
% specifies the result type. '%' for percent values otherswise result is
% W/m^2. Use 'p' for plot and '%p' for percentage plots. The first 5
% parameters are required.
%
% $Id: albedo_altitude.m,v 1.9 2006/05/17 14:39:16 danji Exp $

function result = albedo_altitude(az,pa,alt1,alt2,refl,n,type);

CONST.EMR = 6371.01e3;
CONST.AU = 149598000000;
CONST.AM0 = 1366.9;

if nargin < 5
  error('Not enough input parameters');
end

if nargin < 6
  type = '';
	n = 10;
end

if nargin < 7
  if ischar(n)
    type = n;
    n = 10;
  else
    type = '';
  end
end

maxalt = max(alt1,alt2);
minalt = min(alt1,alt2);
step = floor((maxalt-minalt)/n);
result = zeros(n+1,1);
sat = zeros(3,1);
sun = zeros(3,1);

[sun(1),sun(2),sun(3)] = sph2cart(az,pi/2-pa,CONST.AU);

i = 1;
h = waitbar(0,'Calculating Albedos...');

for alt = minalt:step:maxalt
  [sat(1),sat(2),sat(3)] = sph2cart(az,pi/2-pa,CONST.EMR+alt);
	a = albedo(sat,sun,refl);
	result(i) = sum(sum(a));
	waitbar(i/(n+1),h);
	i = i + 1;
	drawnow;
end

close(h);

if ~isempty(findstr('%',type))
	result = result./CONST.AM0*100;
	yname = 'Irradiance [%]';
else
	yname = 'Irradiance [W/m^2]';
end

if ~isempty(findstr('p',type))
  plot([minalt:step:maxalt]./1000,result);
  title('Subsolar Albedo');
  xlabel('Altitude [km]');
  ylabel(yname);
end

return
