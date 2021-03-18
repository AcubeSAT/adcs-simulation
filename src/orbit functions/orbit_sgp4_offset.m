% -----------------------------------------------------------------------------
%
%                              procedure msrCookieFinal
%
%  This function propagates the satellite by utilizing 
%  the SPG4 orbit propagator software. Based on the propagator, 
%  the function returns magnetic field, sun position and eclipse calculations.
%
%   inputs        :
%     satrec      - struct including all required SPG4 orbit propagator variables
%     dt          - timestep for the orbit propagator
%     ntps        - total simulation seconds
%     tsince_offset - offset time for the orbit propagator (sec)
%
%   outputs       :
%     xsat_ecf      - satellite position in the ECEF frame
%     vsat_ecf      - satellite velocity in the ECEF frame
%     xsat_eci      - satellite position in the ECI frame
%     vsat_eci      - satellite velocity in the ECI frame
%     sat_llh       - satellite position in Latitude, Longitude, Altitude
%     eclipse       - existence or not of eclipse conditions
%     mag_field_ned - reference magnetic field vector in NED frame
%     mag_field_eci - reference magnetic field vector in ECI frame
%     mag_field_ecef - reference magnetic field vector in ECEF frame
%     mag_field_orbit - reference magnetic field vector in Orbit frame
%     sun_pos_ned   - reference sun position vector in NED frame
%     sun_pos_eci   - reference sun position vector in ECI frame 
%     sun_pos_ecef  - reference sun position vector in ECEF frame
%     sun_pos_orbit - reference sun position vector in Orbit frame
%     satrec        - struct including all required SPG4 orbit propagator variables
%     argpm         - argument of perigee 
%     nodem         - right ascension of the ascending node
%     inclm         - inclination                         
%     mm            - mean anomaly                                
%     xnode         - right ascension of the ascending node for short period periodics                 
%     xinc          - inclination for short period periodics     

%  ----------------------------------------------------------------------------*/


function [xsat_ecf, vsat_ecf,xsat_eci,vsat_eci, sat_llh,eclipse, mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit, sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,satrec,argpm,nodem,inclm,mm,xnode,xinc] = orbit_sgp4_offset(satrec,dt,npts,tsince_offset)

%% Parameter initialization
dt=dt/60;
npts = npts/(dt*60);
npts=uint32(npts);

%% Initialize and calculate time
tsince_offset = tsince_offset/60;


t_array=double([0:npts-1]);
tsince=tsince_offset+t_array*dt;    %minutes

Eyear= satrec.epochyr + 2000;
[Emon,Eday,Ehr,Emin,Esec] = days2mdh(Eyear,satrec.epochdays);
UTsec=Ehr*3600+Emin*60+Esec;
% fprintf('EPOCH: YEAR MO  DAY UTSEC \n')
% fprintf('      %5i %2i %4i %5.2f \n\n',Eyear,Emon,Eday,UTsec);

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

%% Calculates position and velocity of the satellite in ECI and ECEF every dt
n=1;

while n<=npts
    [satrec, xsat_eci(:,n), vsat_eci(:,n),argpm(:,n),nodem(:,n),inclm(:,n),mm(:,n),xnode(:,n),xinc(:,n)] = sgp4_2(satrec,tsince(n));
    
    %This segment converts eci coordinates to ecf
    %Compute Greenwich Apparent Siderial Time
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

[sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,mag_field_ned, mag_field_eci,mag_field_ecef,mag_field_orbit] = reference_vectors_calc(sat_llh, time,nodem,inclm,argpm,mm,time_gregorian);

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
end
