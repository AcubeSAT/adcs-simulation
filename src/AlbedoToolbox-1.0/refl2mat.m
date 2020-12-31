% REFL2MAT Convert TOMS reflectivity data to MAT file.
%
% refl2mat(filename [, outdir])
%
% Filename is a string of either a single file, wildcard name, or
% directory. Output files are in outdir or same path as input if not
% specified. Only '.txt' files are read. Data should be loaded with 'refl =
% load(filename)', which will load data into a REFL struct. The output
% filename is of the form gaYYMMDD.mat for compatibility with the
% ALBEDO_WRAPPER simulink interface script.
%
% $Id: refl2mat.m,v 1.12 2006/05/17 14:39:17 danji Exp $

function refl2mat(file, outdir);

if isdir(file)
  files = dir(fullfile(file,'*.txt'));
else
  files = dir(file);
end

if nargin < 2
  if isdir(file)
    outdir = file;
  else
    outdir = fileparts(file);
  end
end

n = size(files,1);

if n == 0
  error('File not found.');
end

for i=1:n
  fprintf('\rFile %i of %i: ',i,n);
  [pathstr,name,ext] = fileparts(files(i).name);
  if strcmp(ext,'.txt')
    refl = read_refl(files(i).name);
    data = refl.data;
    start_time = refl.start_time;
    stop_time = refl.stop_time;
    type = refl.type;
    [startY, startM, startD] = jd2d(refl.start_time);
    if startY-2000 < 0
      startY = startY-1900;
    else
      startY = startY-2000;
    end
    filename = sprintf('ga%02d%02d%02d.mat',startY,startM,startD);
    save(filename,'data','start_time','stop_time','type');
  end
end

fprintf('\n');

return
