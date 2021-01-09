% JD2FILE Returns the filename of the latest REFL data according to a Julian
% Date. The reflectivity data files must be in the path, and have the form
% 'gaYYMMDD.mat' and load a REFL struct named refl. refllib is an optional
% path to the REFL data files, which must contain a directory for each year
% of data. If the data of the specific Julian Date is non-existent, the
% latest file in the same year is returned. If not found the closest future
% date within the year is returned. If still no file is found, the latest
% file from the first previous year is returned, else the first file from
% the closest future year is returned.
%
% filename = jd2file(jd,refllib)
%
% jd is the Julian Date.
%
% $Id: jd2file.m,v 1.4 2006/05/17 11:13:24 danji Exp $

function filename = jd2file(jd,refllib);

[year,month,date] = jd2d(jd);

% Convert year to 2 digit
fyear = year - 2000;
if fyear < 0
  fyear = year - 1900;
end

thefile = sprintf('ga%02.0f%02.0f%02.0f.mat',fyear,month,date);

if isempty(dir(fullfile(refllib,num2str(year),thefile)))
  file = '';
  % Find files in current year
  wfile = sprintf('ga%02.0f*.mat',fyear);
  D = dir(fullfile(refllib,num2str(year),wfile));
  for i=1:size(D,1)
    if sscanf(D(i).name,'ga%f.mat') > sscanf(thefile,'ga%f.mat')
      file = fullfile(refllib,num2str(year),D(i-1).name);
      break;
    end
  end
  % No previous file found. Use the first future file.
  if ~isempty(D) && isempty(file)
    % Filter average etc. files
    for i=1:size(D,1)
      if length(D(i).name) == 12
        file = fullfile(refllib,num2str(year),D(i).name);
        break;
      end
    end
  end
else
  file = fullfile(refllib,num2str(year),thefile);
end

if isempty(file)
  % Find data from previous or future years
  dirs = [];
  D = dir(refllib);
  index = 1;
  firstindex = 0;
  for i=1:size(D,1)
    if D(i).isdir
      dname = str2num(D(i).name); 
      if ~isempty(dname)
        dirs(index) = dname;
        % Mark first index for search
        if dirs(index) < year
          firstindex = index;
        end
        index = index + 1;
      end
    end
  end
  % Search previous years
  for i=firstindex:-1:1
    D = dir(fullfile(refllib,num2str(dirs(i)),'ga*.mat'));
    if ~isempty(D)
      % Filter average etc. files
      for j=size(D,1):-1:1
        if length(D(j).name) == 12
          file = fullfile(refllib,num2str(dirs(i)),D(j).name);
          break;
        end
      end
      % If valid file found then break.
      if ~isempty(file)
        break;
      end
    end
  end
  % Search future years
  for i=firstindex+1:size(dirs,2)
    D = dir(fullfile(refllib,num2str(dirs(i)),'ga*.mat'));
    if ~isempty(D)
      % Filter average etc. files
      for j=1:size(D,1)
        if length(D(j).name) == 12
          file = fullfile(refllib,num2str(dirs(i)),D(j).name);
          break;
        end
      end
    end
    % If valid file found then break.
    if ~isempty(file)
      break;
    end
  end
end

if isempty(file)
  error('jd2file: No reflectivity data found.');
else
  filename = file;
end

return
