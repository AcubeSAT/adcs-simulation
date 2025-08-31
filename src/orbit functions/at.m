
filename = 'new_thomson.a'; 

% This script generates an attitude file for use with STK.
%It is to be executed only after the profile simulation.

% Inputs(modify as needed based on the selected profile):
% Time                - N x 1 matrix (in seconds)
% quaternion          - 4 x N matrix 
% Angular Velocities  - 3 x N matrix  (in deg/sec)

numPoints = length(Time);  % N = number of points


% Load TLE from a text file (adjust file name/path as needed)
tle_file = "\TLEs_2027\SSO-600-11PM.TLE"; 
fileID = fopen(tle_file, 'r');
if fileID == -1
    error('Cannot open TLE file: %s', tle_file);
end

% Read the two lines of TLE
tle_lines = cell(1,2);
tle_lines{1} = fgetl(fileID);
tle_lines{2} = fgetl(fileID);
fclose(fileID);


% Parse epoch from TLE line 1 (YYDDD.DDDDDDDD format)
tleEpochStr = tle_lines{1}(19:32);
year = str2double(tleEpochStr(1:2));
doy = str2double(tleEpochStr(3:end));

if year < 57
    year = year + 2000;
else
    year = year + 1900;
end


epochDate = datetime(year, 1, 1) + days(doy - 1);
epochStr = datestr(epochDate, 'dd mmm yyyy HH:MM:SS');  


% ========== Write STK .a File ========== 
fid = fopen(filename, 'w');

fprintf(fid, 'stk.v.12.2\n\n');
fprintf(fid, 'BEGIN Attitude\n\n');
fprintf(fid, '    NumberOfAttitudePoints %d\n', numPoints);
fprintf(fid, '    CentralBody   Earth\n');
fprintf(fid, '    ScenarioEpoch   %s\n\n', epochStr);
fprintf(fid, '    CoordinateAxes ICRF\n\n');
fprintf(fid, 'AttitudeTimeQuatAngVels\n\n');

rad_to_deg = 180 / pi; %conversion from  radians to degrees 


ang_vel=zeros(3,length(Time));
q=zeros(4,length(Time));
for i = 1:numPoints
    % The format must be as follows:
    % [Time Quat(1) Quat(2) Quat(3) Quat(4) AngVelX AngVelY AngVelZ]
R=[0 -1 0;0 0 -1;1 0 0];
disp(R);
ang_vel(:,i)=R*x_real(5:7,i);
quat_rot=d2q(R);
q(:,i)=quatProd(quat_rot,x_real(:,i));
%ang_vel(:,i)=x_real(5:7,i)'*R;
    fprintf(fid, '%.3f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
        Time(i), q(1,i), q(2, i), q(3, i), q(4, i), ...
      ang_vel(1,i)*rad_to_deg,ang_vel(2,i)*rad_to_deg, ang_vel(3,i)*rad_to_deg);
end

fprintf(fid, '\nEND Attitude\n');
fclose(fid);

fprintf('STK attitude file written to "%s"\n', filename);
