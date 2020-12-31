% LOAD_REFL Reads REFL MAT file and automatically scans the refl_data
% directory, after current dir, in the Albedo Toolbox installation
% directory.
%
% refl = load_refl(date_str)
%
% where date is a text string specifying the date of the data in the yymmdd
% format or yymmdd-yymmdd format for statistical calculations. Current dir
% is searched first, but only if date_str contains the full filename with
% or without .mat extension.
%
% $Id: load_refl.m,v 1.6 2006/05/17 14:39:17 danji Exp $

function refl = load_refl(date_str)

%%%%%%%%%%%%%%%%%%
% Init
%%%%%%%%%%%%%%%%%%

if nargin < 1
    % No date specified
    error('No date specified.');
end

filename = strcat('ga',date_str,'.mat');

%%%%%%%%%%%%%%%%%%
% Try current dir
%%%%%%%%%%%%%%%%%%

% Try if full filename was specified

if size(dir(date_str),1) > 0
  refl = load(date_str);
  return
end

% Try if filename without .mat extension was specified

if size(dir(strcat(date_str,'.mat')),1) > 0
  refl = load(strcat(date_str,'.mat'));
  return
end

%%%%%%%%%%%%%%%%%%
% Find year
%%%%%%%%%%%%%%%%%%

YY  = date_str(1:2);

if double(YY) > 70
  year = strcat('19',YY);
else
  year = strcat('20',YY);
end

%%%%%%%%%%%%%%%%%%
% Load file
%%%%%%%%%%%%%%%%%%

file = fullfile(albedo_path,'refl_data',year,filename);

try
  refl = load(file);
catch
  lasterr_struct = lasterror;
  msg = sprintf('File not found for date %s.\n\n??? %s',date_str,lasterr_struct.message);
  error(msg);
end

return
