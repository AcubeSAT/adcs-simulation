% JD2D Convert Julian Date to Gregorian Calendar Date.
%
% [year,month,day,hour,min,sec] = jd2d(jd)
%
% Precision can be reduced by leaving out output parameters from right to
% left, e.g. [year,month] = jd2d(jd) returns the year and month of jd.
% Resolution is in milliseconds.
%
% $Id: jd2d.m,v 1.2 2006/05/17 14:39:18 danji Exp $

function [year,month,day,hour,min,sec] = jd2d(jd);

Z = floor(jd+0.5);
F = jd+0.5-Z;

if  Z < 2299161
  A = Z;
else
  alpha = floor((Z-1867216.25)/36524.25);
  A = Z+1+alpha-floor(alpha/4);
end

B = A+1524;
C = floor((B-122.1)/365.25);
D = floor(365.25*C);
E = floor((B-D)/30.6001);

dd = B-D-floor(30.6001*E)+F;
if E < 13.5
  mm = E-1;
else
  mm = E-13;
end
if mm > 2.5
  year = C-4716;
else
  year = C-4715;
end

if nargout > 1
  month = mm;
end
if nargout > 2
  day = floor(dd);
end
if nargout > 3
  hh = 24*(dd-day);
  hour = floor(hh);
end
if nargout > 4
  ii = 60*(hh-hour);
  min = floor(ii);
end
if nargout > 5
  sec = 60*(ii-min);
  % Scale to msec precision
  sec = round(sec*1000)/1000;
end

return
