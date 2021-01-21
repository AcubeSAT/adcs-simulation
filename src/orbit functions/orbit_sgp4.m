function [xsat_ecf, vsat_ecf,xsat_eci,vsat_eci, sat_llh,eclipse, mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit, sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,satrec,argpm,nodem,inclm,mm,xnode,xinc] = orbit_sgp4(satrec,dt,npts)

% dt is xronos deigmatolipsias
% npts is times the loop must be run (orbits*orbitPeriod)

%% Initializations
dt=dt/60;
npts = npts/(dt*60);
npts=uint32(npts);

%% Initialize and calculate time
tsince_offset =2000/60;

t_array=double([0:npts-1]);
tsince=tsince_offset+t_array*dt;%minutes

Eyear= satrec.epochyr + 2000;
[Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,satrec.epochdays);
UTsec=Ehr*3600+Emin*60+Esec;
fprintf('EPOCH: YEAR MO  DAY UTSEC \n')
fprintf('      %5i %2i %4i %5.2f \n\n',Eyear,Emon,Eday,UTsec);

time_gregorian=zeros(length(tsince),1);
time_day=satrec.epochdays+tsince/1440;

for i=1:length(tsince)
    if time_day(i)>=366
        Eyear=Eyear+1;
        time_day=time_day-365;
    end
    [mon,day,hr,min,sec] = days2mdh(Eyear,time_day(i));
    greg_date_vector=[Eyear, mon, day, hr, min, sec];

    time_gregorian(i)=datenum(greg_date_vector);
end
time_gregorian=time_gregorian';

%% Calculates x and v of the satellite in ECI and ECEF every dt
n=1;

while n<=npts
    [satrec, xsat_eci(:,n), vsat_eci(:,n),argpm(:,n),nodem(:,n),inclm(:,n),mm(:,n),xnode(:,n),xinc(:,n)] = sgp4_2(satrec,tsince(n));
    
    %This segment converts eci coordinates to ecf
    %Compute Greenwich Apparent Siderial Time
    tsince(n);
    gst=gstime(satrec.jdsatepoch+tsince(n)/1440);
    [xsat_ecf(:,n),vsat_ecf(:,n)] = ECI2ECEF(xsat_eci(:,n),vsat_eci(:,n),gst,1);
    n=n+1;
end

xsat_ecf=xsat_ecf*1000;  %m
vsat_ecf=vsat_ecf*1000;  %mps
xsat_eci=xsat_eci*1000;  %m
vsat_eci=vsat_eci*1000;  %mps

%% Calculate llh (rad)
sat_llh =ecf2llhT(xsat_ecf);

%% Calculate reference vectors
time = satrec.jdsatepoch + tsince(1:npts)/1440;

[sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,mag_field_ned, mag_field_eci,mag_field_ecef,mag_field_orbit] = reference_vectors_calc(sat_llh, time,nodem,inclm,argpm,mm,tsince,time_gregorian);

%% Calculate Eclipse
eclipse = calculate_eclipse(xsat_eci/1000,sun_pos_eci);

penumbral=0;
umbral=0;
normal=0;
for i=1:size(eclipse,2)
if eclipse(i) == 1
    penumbral=penumbral+1;
elseif eclipse(i) == 2
    umbral=umbral+1;
else
    normal=normal+1;
end
end
normal;
penumbral;
umbral;
normal+penumbral+umbral;
npts;


end