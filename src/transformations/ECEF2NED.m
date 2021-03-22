% -----------------------------------------------------------------------------
%
%                              procedure ECEF2NED
%
%  this procedure rotates a vector from the ECEF to the NED frame
%
%   inputs        :
%   ref_vec_ecef  - Vector in ECEF frame
%   lat           - Latitude
%   long          - Longitude
%
%   outputs       :
%   ref_vector_ned - Vector in NED frame
%  ----------------------------------------------------------------------------*/

function ref_vector_ned = ECEF2NED(ref_vec_ecef,lat,long)
for n=1:size(ref_vec_ecef,2)
    c1=[-sin(lat(n))*cos(long(n));-sin(lat(n))*sin(long(n));cos(lat(n))];
    c2=[-sin(long(n));cos(long(n));0];
    c3=[-cos(lat(n))*cos(long(n));-cos(lat(n))*sin(long(n));-sin(lat(n))];
    R = [c1,c2,c3];
    ref_vector_ned(:,n)=R'*ref_vec_ecef(:,n);
end
end
