% ALBEDO_WRAPPER Simulink wrapper for ALBEDO.
%
% a = albedo_wrapper(sat, sun, param, redfac [, refllib])
%
% sat and sun are the vectors to the Earth and Sun in Earth Centered Earth
% Fixed coordinates, respectively. Param is either a REFL struct or a
% julian date. redfac is the reduction factor, which is used to reduce the
% resolution of the REFL data. refllib is the path to the reflectivity data,
% only needed if param is a Julian Date.
%
% $Id: albedo_wrapper.m,v 1.14 2006/05/17 14:39:16 danji Exp $

function a = albedo_wrapper(sat,sun,param,redfac,refllib);

CONST.AM0 = 1366.9;
persistent refl;
persistent rrefl;
persistent oldredfac;

% Reload refl if reduction factor is changed
if redfac ~= oldredfac
  refl = [];
end
oldredfac = redfac;

if isstruct(param)
  if ~isstruct(rrefl) || rrefl.start_time ~= param.start_time || rrefl.stop_time ~= param.stop_time
    rrefl = resize_refl(param,redfac);
  end
  a = albedo(sat,sun,rrefl);
else
  if ~(isfield(refl,'data') && param >= refl.start_time && param <= refl.stop_time)
    refl = load(jd2file(param,refllib));
    refl = replace_nan(refl,refllib);
    refl = resize_refl(refl,redfac);
    % Check if data is unavailable. If so, earliest data is selected, and
    % is valid from current time to stop time.
    if param < refl.start_time
      refl.start_time = param;
    end
  end
  a = albedo(sat,sun,refl);
end

return
