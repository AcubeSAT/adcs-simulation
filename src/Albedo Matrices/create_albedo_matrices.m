%% Script for creating albedo matrices
%
% Input: 
% This script requires the AlbedoToolbox-1.0 and the refl_data folder with
% the reflectivity data for a given date.
%
% Output: 
% Albedo as a percentage of sun irradiance for the duration defined by running_time
% with timestep dt.
%
% Note: NaN values are replaced by the annual average.

reflectivity = load_refl('050403');
dt = 1;
running_time =(5545+0.9);
satrec = orbit_init();

[xsat_ecf, vsat_ecf,xsat_eci,vsat_eci, sat_llh,eclipse, mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit, sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,satrec,argpm,nodem,inclm,mm,xnode,xinc] = orbit_sgp4(satrec,dt,running_time+dt);

new_refl = replace_nan(reflectivity,'ga050101-051231.mat');

sun_pos_ecef_m = 1.495978707e+11 * sun_pos_ecef;

total_albedo = zeros(1,length(xsat_ecf));

for i = 1:length(xsat_ecf)
    a = albedo(xsat_ecf(:,i),sun_pos_ecef_m(:,i),new_refl);
    total_albedo(i) = sum(sum(a)); 
end

fine_grained_total_albedo = zeros(1,length(total_albedo)*10);
for i=1:length(total_albedo)
    fine_grained_total_albedo((i-1)*10+1:i*10 ) = total_albedo(i);
end

