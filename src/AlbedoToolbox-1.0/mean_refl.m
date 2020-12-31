% MEAN_REFL Calculate mean of REFL data in MAT files.
%
% [mean_val, mean_lat] = mean_refl(files)
%
% files is a string of a wildcard filename or a directory. Only files with
% '.mat' extension are read. mean_val contains the mean values and mean_lat
% contains the mean values meaned over latitudes. mean_val.start_time is the
% earliest file start time, and mean_val.stop_time is the latest. mean_val
% is a REFL struct. mean_lat is a vector.
%
% $Id: mean_refl.m,v 1.12 2006/05/17 11:13:23 danji Exp $

function [mean_val, mean_lat] = mean_refl(file)

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

% Allocate valid sample count
refl_acc = zeros(180,288);

% Allocate mean value array
refl_count = zeros(180,288);

% Initialize start and stop times
start_time = zeros(nfiles,1);
stop_time = zeros(nfiles,1);

%%%%%%%%%%%%%%%%%%
% Calculate mean
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
    % Get REFL data
    for index = 1:180*288
      if ~isnan(refl.data(index))
        % Accumelate data values if valid
        refl_acc(index) = refl_acc(index) + refl.data(index);
        % Track valid sample count for mean value calculation
        refl_count(index) = refl_count(index) + 1;
      end
    end
  else
    start_time(ifile) = inf;
  end
  fprintf(1,'\b\b\b\b%3.0f%%',ifile*100/nfiles);
end

% Supress Divide by Zero warnings in case of missing data points 
warning('off');

 % Return mean_val
mean_val = refl_struct(refl_acc./refl_count,min(start_time),max(stop_time),'Mean');
if min(min(refl_count)) == 0
  fprintf(1,'\nmean_val contains %.0f%% empty values (NaN).',sum(sum(isnan(mean_val.data)))*100/(180*288));
end

% Return mean_lat
if nargout > 1
  mean_lat = sum(refl_acc')'./sum(refl_count')';
  if min(sum(refl_count')) == 0
    fprintf(1,'\nmean_lat contains %.0f%% empty values (NaN).',sum(isnan(mean_lat))*100/180);
  end
end

fprintf('\n');

return
