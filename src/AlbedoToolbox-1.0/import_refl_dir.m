% IMPORT_REFL_DIR Imports TOMS reflectivity data to the Earth Albedo
% Toolbox. REFL files with suffix '.txt' are converted using
% refl2mat, and the statistics are calculated using MEAN_REFL and STD_REFL,
% and saved in the directory.
%
% import_refl_dir(refldir)
%
% refldir is a string containing the path to the directory holding the
% REFL files.
%
% $Id: import_refl_dir.m,v 1.5 2006/05/17 14:39:17 danji Exp $

function import_refl_dir(refldir)

matfiles = dir(fullfile(refldir,'*.mat'));

if size(matfiles,1) > 0
  ustr = '';
  while ~strcmp(lower(ustr),'n') && ~strcmp(lower(ustr),'y')
    ustr = input('\nMat-files exist. Convert .txt files? [y/n]: ', 's');
  end
  if strcmp(lower(ustr),'y')
    fprintf(1, 'Converting REFL files:\n');
    refl2mat(refldir);
  end
else
    fprintf(1, 'Converting REFL files:\n');
    refl2mat(refldir);
end

fprintf(1, 'Calculating mean:\n');
[mean_val, mean_lat] = mean_refl(refldir);
fprintf(1, 'Calculating standard deviation:\n');
[std_val, std_lat] = std_refl(refldir,mean_val,mean_lat);

filename = ['ga' get_date_str(mean_val) '.mat'];
filename_lat = ['ga' get_date_str(mean_val) 'lat.mat'];
data = mean_val.data;
start_time = mean_val.start_time;
stop_time = mean_val.stop_time;
type = mean_val.type;
save(filename, 'data', 'start_time', 'stop_time', 'type');
fprintf(1,'Saved %s\n',filename);
save(filename_lat, 'mean_lat');
fprintf(1,'Saved %s\n',filename_lat);

filename = ['ga' get_date_str(mean_val) 'std.mat'];
filename_lat = ['ga' get_date_str(mean_val) 'latstd.mat'];
data = std_val.data;
start_time = std_val.start_time;
stop_time = std_val.stop_time;
type = std_val.type;
save(filename, 'data', 'start_time', 'stop_time', 'type');
fprintf(1,'Saved %s\n',filename);
save(filename_lat, 'std_lat');
fprintf(1,'Saved %s\n',filename_lat);

fprintf(1,'Complete. Original REFL files are not deleted.\n');

return

function datestr = get_date_str(refl);
% Construct datestamp for filenames

[startY, startM, startD] = jd2d(refl.start_time);
% Take the last whole day (prevent stop time to be 0:00 on the day after)
[stopY, stopM, stopD] = jd2d(refl.stop_time-1/24);

% Y2K

if startY-2000 < 0
  startY = startY-1900;
else
  startY = startY-2000;
end  
if stopY-2000 < 0
  stopY = stopY-1900;
else
  stopY = stopY-2000;
end  

startstr = sprintf('%02d%02d%02d',startY,startM,startD);
stopstr = sprintf('%02d%02d%02d',stopY,stopM,stopD);

datestr = [startstr '-' stopstr];

return
