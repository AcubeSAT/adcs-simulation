% READ_REFL Reads EP/TOMS REFL data from file into an array. The return value
% is a REFL struct. See 'help refl_struct' for details. 
%
% refl = read_refl(reflfile)
%
% where filename is a text string specifying the file to read.
%
% $Id: read_refl.m,v 1.7 2006/05/17 14:39:17 danji Exp $

function refl = read_refl(reflfile)

%%%%%%%%%%%%%%%%%%
% Init
%%%%%%%%%%%%%%%%%%

if nargin < 1
    % No file specified
    error('No file specified.');
end

[fid msg] = fopen(reflfile,'r');

if fid == -1
  error(msg);
end

%%%%%%%%%%%%%%%%%%
% Read header
%%%%%%%%%%%%%%%%%%

% Read date line
start_time = get_date(fid);
stop_time = start_time+1;

% Discard header line 2 and 3
junk = fgetl(fid);
junk = fgetl(fid);

%%%%%%%%%%%%%%%%%%
% Read data
%%%%%%%%%%%%%%%%%%

data = zeros(180,288);
fprintf(1,'Progress: %3.0f%%',0);

% Latitude from -89.5 to 89.5 deg in 1 deg steps
for ilat = 1:180
  % Read empty space at beginning of every longitude block
  junk = fscanf(fid,'%c',1);
    % Longitude from -179.375 to 179.375 deg in 1.25 deg steps
    for ilon = 1:288
        % Read three digits
        val_str = fscanf(fid,'%c',3);
        % Convert to fraction from percentage characters
        val = str2double(val_str)*0.01;
        % Save data in output matrix
        if val == 9.99
          data(ilat,ilon) = NaN;
        else
          data(ilat,ilon) = val;
        end
        % Read LF [+CR] and space at BOL for every 25 values
        if rem(ilon,25) == 0
          junk = fgetl(fid);
          junk = fscanf(fid,'%c',1);
        end
        % Loop for next value
    end
    % Read block end
    junk = fgetl(fid);
    % Print status
    fprintf(1,'\b\b\b\b%3.0f%%',ilat*100/180);
end

fclose(fid);

refl = refl_struct(data,start_time,stop_time,'Raw');

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Julian Date from header line 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jd = get_date(fid);

% Read Day number
junk = fscanf(fid,'%c',10);

% Read Month
switch fscanf(fid,'%c',3);
  case 'Jan'
    M = 1;
  case 'Feb'
    M = 2;
  case 'Mar'
    M = 3;
  case 'Apr'
    M = 4;
  case 'May'
    M = 5;
  case 'Jun'
    M = 6;
  case 'Jul'
    M = 7;
  case 'Aug'
    M = 8;
  case 'Sep'
    M = 9;
  case 'Oct'
    M = 10;
  case 'Nov'
    M = 11;
  case 'Dec'
    M = 12;
  otherwise
    error('Error reading month in refl file');
end

junk = fscanf(fid,'%c',1);

% Read day
d1 = fscanf(fid,'%c',1);
if d1 == ' '
    d1 = '0';
end

d2 = fscanf(fid,'%c',1);

D = str2num([d1 d2]);

junk = fscanf(fid,'%c',2);

% Read year
Y = str2num(fscanf(fid,'%c',4));

% Read rest of line
junk = fgetl(fid);

% Calculate Julian Date (Cite: Wolfram Astronomy)
jd = d2jd(Y,M,D);

return
