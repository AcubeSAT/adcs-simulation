% -----------------------------------------------------------------------------
%
%                              procedure reference_vectors_calc
%
%  This function calculates the magnetic field and sun position
%  in every ADCS defined frame, for every orbit propagator time and position pair.
%
%   inputs        :
%     sat_llh     - satellite position in Latitude, Longitude, Altitude
%     time        - time in JD format
%     nodeo       - right ascension of the ascending node
%     inclo       - inclination  
%     argpo       - argument of perigee                        
%     mo          - mean anomaly  
%     time_gregorian -  time in gregorian format
%
%   outputs       :
%     mag_field_ned - reference magnetic field vector in NED frame
%     mag_field_eci - reference magnetic field vector in ECI frame
%     mag_field_ecef - reference magnetic field vector in ECEF frame
%     mag_field_orbit - reference magnetic field vector in Orbit frame
%     sun_pos_ned   - reference sun position vector in NED frame
%     sun_pos_eci   - reference sun position vector in ECI frame 
%     sun_pos_ecef  - reference sun position vector in ECEF frame
%     sun_pos_orbit - reference sun position vector in Orbit frame
%
%  ----------------------------------------------------------------------------*/


function [sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit] = reference_vectors_calc(sat_llh, time,nodeo,inclo,argpo,mo,time_gregorian)

%% Calculate Magnetic reference field

% mag_field_ned=zeros(3,length(time));
% for n=1:size(time,2)
%     mag_field_ned(:,n) = igrf(time_gregorian(n),sat_llh(1,n)*180/pi,sat_llh(2,n)*180/pi,sat_llh(3,n)/1000,'geodetic');
% end

% uncomment the three following lines for speed purposes (sets the time constant - minor calculation mistake)
time_init_greg=time_gregorian(1);
mag_field_ned = igrf(time_init_greg,sat_llh(1,:)*180/pi,sat_llh(2,:)*180/pi,sat_llh(3,:)/1000,'geodetic');
mag_field_ned=mag_field_ned';

%% Find the Magnetic Field in all used frame
mag_field_ecef = NED2ECEF(mag_field_ned,sat_llh(1,:),sat_llh(2,:));
mag_field_eci=zeros(3,length(time));
mag_field_orbit=zeros(3,length(time));
for n=1:size(mag_field_ecef,2)
    
    gst=gstime(time(n));
    % gst1=JD2GAST(time(n))*pi/180;
    [mag_field_eci(:,n),R] = ECEF2ECI(mag_field_ecef(:,n),zeros(3,1),gst,0);
    mag_field_orbit(:,n) = ECI2Orbit(mag_field_eci(:,n),nodeo(1,n),inclo(1,n),argpo(1,n)+mo(1,n));

end


%% Sun position reference vector
sun_pos_eci=sun_position(time);
sun_pos_orbit=zeros(3,length(time));
sun_pos_ecef=zeros(3,length(time));
for n=1:size(mag_field_eci,2)
    sun_pos_orbit(:,n)=ECI2Orbit(sun_pos_eci(:,n),nodeo(1,n),inclo(1,n),argpo(1,n)+mo(1,n));
    gst=gstime(time(n));
    [sun_pos_ecef(:,n),] = ECI2ECEF(sun_pos_eci(:,n),zeros(3,1),gst,0);
end
sun_pos_ned=ECEF2NED(sun_pos_ecef,sat_llh(1,:),sat_llh(2,:));

end
