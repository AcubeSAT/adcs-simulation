% STD_REFL Calculate standard deviation of REFL data in MAT files.
%
% [std_val, std_lat] = std_refl(files, mean_val, mean_lat)
%
% files is a string of a wildcard filename or a directory. Only files with
% '.mat' extension are read. mean_val and mean_lat are the mean values of
% the REFL data in files. std_val contains the standard deviation values
% and std_lat contains the standard deviation values over latitudes.
% std_val.start_time is the earliest file start time, and std_val.stop_time
% is the latest. std_val is a REFL struct. std_lat is a vector.
%
% $Id: std_refl.m,v 1.15 2006/05/17 14:39:17 danji Exp $

function [std_val, std_lat] = std_refl(file, mean_val, mean_lat)

%%%%%%%%%%%%%%%%%%
% Check for existing statistics files
%%%%%%%%%%%%%%%%%%

statfiles = dir('*-*.mat');

if size(statfiles,1) > 0
    error('Please remove mat files containing mean or standard deviation, so they are not included in the calculations');
    return
end

%%%%%%%%%%%%%%%%%%
% Init
%%%%%%%%%%%%%%%%%%

[parm_path,name,ext,versn] = fileparts(file);

if isdir(file)
  parm_path = file;
  files = dir(fullfile(file,'*.mat'));
else
  files = dir(file);
end

nfiles = size(files,1);

% Allocate accumulated variance array
refl_acc = zeros(180,288);

% Allocate accumulated variance array with 180x1 mean array for simplified
% latitude based refl model
refl_lat_acc = zeros(180,288);

% Allocate valid sample count
refl_count = zeros(180,288);

% Initialize start and stop times
start_time = zeros(nfiles,1);
stop_time = zeros(nfiles,1);

%%%%%%%%%%%%%%%%%%
% Calculate standard deviation
%%%%%%%%%%%%%%%%%%

fprintf(1,'Progress: %3.0f%%',0);

for ifile = 1:nfiles
  [pathstr,name,ext,vers] = fileparts(files(ifile).name);
  if strcmp(ext,'.mat')
    % Read next file
    refl = load(fullfile(parm_path, files(ifile).name));
    % Get start and stop times
    start_time(ifile) = refl.start_time;
    stop_time(ifile) = refl.stop_time;
    for index = 1:180*288
      if ~isnan(refl.data(index))
        % Accumelate data values if valid
        refl_acc(index) = refl_acc(index) + (refl.data(index)-mean_val.data(index))^2;
        refl_lat_acc(index) = refl_lat_acc(index) + (refl.data(index)-mean_lat(mod(index-1,180)+1))^2;
        % Track valid sample count for mean value calculation
        refl_count(index) = refl_count(index) + 1;
      end
    end
    start_time(ifile) = inf;
  end
  fprintf(1,'\b\b\b\b%3.0f%%',ifile*100/nfiles);
end

% Return std_val
warning('off');
std_val = refl_struct(sqrt(refl_acc./refl_count),min(start_time),max(stop_time),'Standard Deviation');
if min(min(refl_count)) == 0
  fprintf(1,'\nstd_val contains %.0f%% empty values (NaN).',sum(sum(isnan(std_val.data)))*100/(180*288));
end

% Return std_lat
if nargout > 1
  std_lat = sqrt(sum(refl_lat_acc')./sum(refl_count'))';
  if min(sum(refl_count')) == 0
    fprintf(1,'\nstd_lat contains %.0f%% empty values (NaN).',sum(isnan(std_lat))*100/180);
  end
end

fprintf('\n');

return
