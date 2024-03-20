% -----------------------------------------------------------------------------
%
%                              procedure ECEF2NED2
%
%  this procedure rotates a vector from the ECEF to the NED frame
%
%  Inputs          :
%   ref_vec_ecef   - Vector in ECEF frame
%   lat            - Latitude
%   long           - Longitude
%
%  Outputs         :
%   ref_vector_ned - Vector in NED frame
%
%
%  Based on AOCS_DDJF
%  ----------------------------------------------------------------------------*/

function ref_vector_ned = ECEF2NED(ref_vec_ecef,lat,long)

% Construct rotation matrix
c1=[-sin(lat)*cos(long); -sin(lat)*sin(long); cos(lat)];
c2=[-sin(long)         ; cos(long)          ; 0];
c3=[-cos(lat)*cos(long); -cos(lat)*sin(long); -sin(lat)];

% Transpose rotation matrix
R = [c1,c2,c3]';

% Rotate vector
ref_vector_ned = R*ref_vec_ecef;

end
