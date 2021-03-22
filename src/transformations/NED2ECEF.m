% -----------------------------------------------------------------------------
%
%                              procedure NED2ECEF
%
%  this procedure rotates a vector from the NED to the ECEF frame
%
%   inputs        :
%   ref_vec_ned   - Vector in NED frame
%   lat           - Latitude
%   long          - Longitude
%
%   outputs       :
%   ref_vector_ecef - Vector in ECEF frame
%  ----------------------------------------------------------------------------*/
function ref_vector_ecef = NED2ECEF(ref_vec_ned,lat,long)
for n=1:size(ref_vec_ned,2)
    c1=[-sin(lat(n))*cos(long(n));-sin(lat(n))*sin(long(n));cos(lat(n))];
    c2=[-sin(long(n));cos(long(n));0];
    c3=[-cos(lat(n))*cos(long(n));-cos(lat(n))*sin(long(n));-sin(lat(n))];
    R = [c1,c2,c3];
    ref_vector_ecef(:,n)=R*ref_vec_ned(:,n);
end
end
