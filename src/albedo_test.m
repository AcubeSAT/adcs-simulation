refl = load_refl('050403');
dt = 1;
tf =(5545+0.9);
satrec = orbit_init();

[xsat_ecf, vsat_ecf,xsat_eci,vsat_eci, sat_llh,eclipse, mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit, sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,satrec,argpm,nodem,inclm,mm,xnode,xinc] = orbit_sgp4(satrec,dt,tf+dt);

new_refl = replace_nan(refl,'ga050101-051231.mat');

sun_pos_ecef_m = 1.495978707e+11 * sun_pos_ecef;


total_albedo = zeros(1,length(xsat_ecf));

for i = 1:length(xsat_ecf)
    a = albedo(xsat_ecf(:,i),sun_pos_ecef_m(:,i),new_refl);
    total_albedo(i) = sum(sum(a)); % summing all elements
end

new = zeros(1,length(total_albedo)*10);
for i=1:length(total_albedo)
    new((i-1)*10+1:i*10 ) = total_albedo(i);
end

