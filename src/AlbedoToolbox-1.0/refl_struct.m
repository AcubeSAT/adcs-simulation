% REFL_STRUCT Create REFL struct from parameters. The struct contains the following
% fields:
%    - data (reflectivity data)
%    - start_time (Julian Date of start time)
%    - stop_time (Julian Date of stop time)
%    - type (String specifying data description, e.g. Raw or Mean)
%
% refl = refl_struct(data,start_time,stop_time,type)
%
% $Id: refl_struct.m,v 1.2 2006/02/23 08:31:33 danji Exp $

function refl = refl_struct(data,start_time,stop_time,type);

refl.data = data;
refl.start_time = start_time;
refl.stop_time = stop_time;
refl.type = type;

return
