
filename = 'attitude.a'; 

% This script generates an attitude file for use with STK.
%It is to be executed only after the profile simulation.

% Inputs(modify as needed based on the selected profile):
% Time                - N x 1 matrix (in seconds)
% quaternion          - 4 x N matrix 
% Angular Velocities  - 3 x N matrix  (in deg/sec)

numPoints = length(Time);  % N = number of points

% ========== Get Scenario Epoch from TLE ========== 
tle_lines = {'1 99999U          27112.41666667  .00000000  00000-0  00000-0 0 00003', ...
             '2 99999 097.4116 195.2397 0011825 245.9291 114.0710 15.24069409000015'};

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




for i = 1:numPoints
    % The format must be as follows:
    % [Time Quat(1) Quat(2) Quat(3) Quat(4) AngVelX AngVelY AngVelZ]
    
    fprintf(fid, '%.3f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
        Time(i), q_ob_data(1,i), q_ob_data(2, i), q_ob_data(3, i), q_ob_data(4, i), ...
        x_real(5, i)*rad_to_deg, x_real(6, i)*rad_to_deg, x_real(7, i)*rad_to_deg);
end

fprintf(fid, '\nEND Attitude\n');
fclose(fid);

fprintf('STK attitude file written to "%s"\n', filename);
