% D2JD Convert Gregorian Calendar Date to Julian Date.
%
% jd = d2jd(year,month,day,hour,min,sec)
%
% Precision can be reduced by leaving out parameters from right to left,
% e.g. d2jd(2004) gives the Julian Date of Jan 1, 2004 at 0:00 am.
%
% $Id: d2jd.m,v 1.2 2006/05/17 14:39:17 danji Exp $

function jd = d2jd(year,month,day,hour,min,sec);

if nargin < 6
  sec = 0;
end
if nargin < 5
  min = 0;
end
if nargin < 4
  hour = 0;
end
if nargin < 3
  day = 1;
end
if nargin < 2
  month = 1;
end
if nargin < 1
  error('At least one parameter is needed.');
end

jd = 367*year-floor(7*(year+floor((month+9)/12))/4) ...
  -floor(3*(floor((year+(month-9)/7)/100)+1)/4) ...
  +floor(275*month/9)+day+1721028.5+(hour+(min+sec/60)/60)/24;

return
