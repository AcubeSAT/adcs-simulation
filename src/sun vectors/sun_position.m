% -----------------------------------------------------------------------------
%
%                              procedure sun_position
%
%  This procedure estimates the predicted sun position on the
%  Earth-Centered Inertial frame based on the satellite time. 
%
%   inputs        :
%   time          - Time in JD format
%
%   outputs       :
%     sun_pos_eci - Sun position in the ECI frame
%  ----------------------------------------------------------------------------*/

function [sun_pos_eci] = sun_position(time)

% Convert time
ut1=(time - 2451545.0) / 36525.0;

%% ------ Algorithm --------
meanlong = 280.4606184 + 36000.77005361*ut1;
meanlong=mod(meanlong,360); %deg
meananomaly=357.5277233+35999.05034*ut1;
meananomaly=mod(meananomaly*pi/180,2*pi); %rad
if meananomaly < 0
   meananomaly = 2*pi + meananomaly;
end
eclplong = meanlong + 1.91466471*sin(meananomaly) + 0.019994643 * sin(2 * meananomaly); %deg
obliquity = 23.439291 - 0.0130042 * ut1; %deg
meanlong = meanlong*pi/180; %rad
if meanlong < 0
    meanlong = 2*pi + meanlong;
end
eclplong=eclplong*pi/180; %rad
obliquity = obliquity*pi/180; %rad

%% Magnitude of sun vector and ECI calculations
magr = 1.000140612 - 0.016708617*cos(meananomaly)-0.000139589*cos(2*meananomaly); %au's
sun_pos_eci = [magr.*cos(eclplong); magr.*cos(obliquity).*sin(eclplong); magr.*sin(obliquity).*sin(eclplong)];


rtasc = atan(cos(obliquity).*tan(eclplong));
if eclplong < 0 
    eclplong = eclplong +2*pi;
end

%% check rtasc in the same quadrant as eclplong
for n=1:length(eclplong)
if abs(eclplong(n) - rtasc(n)) > (pi/2)
    if mod(eclplong(n) - rtasc(n),1)<0.51
        rtasc(n) = rtasc(n) + 0.5*pi*((eclplong(n) - rtasc(n))-mod(eclplong(n) - rtasc(n),1))/(pi/2);
    else
        rtasc(n) = rtasc(n) + 0.5*pi*((eclplong(n) - rtasc(n))+(1-mod(eclplong(n) - rtasc(n),1)))/(pi/2);
    end
end
end

decl = asin(sin(obliquity).*sin(eclplong));
end
