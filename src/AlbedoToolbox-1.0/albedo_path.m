% ALBEDO_PATH Returns the path to albedo_path.m, under the assumption that this% is the albedo toolbox installation dir.
%
%  pathstr = albedo_path
%
% $Id: albedo_path.m,v 1.3 2006/05/17 14:39:16 danji Exp $

function pathstr = albedo_path;

[pathstr,name,ext] = fileparts(which('albedo_path.m'));

return
