% -----------------------------------------------------------------------------
%
%                              procedure NED2ECEF
%
%  this procedure rotates a vector from the NED to the ECEF frame
%
%  Inputs           :
%   ref_vec_ned     - Vector in NED frame
%   lat             - Latitude
%   long            - Longitude
%
%  Outputs          :
%   ref_vector_ecef - Vector in ECEF frame
%
%
%  Based on AOCS_DDJF
%  ----------------------------------------------------------------------------*/

function ref_vector_ecef = NED2ECEF(ref_vec_ned, lat, long)

    % % Construct rotation matrix
    c1 = [-sin(lat) * cos(long); -sin(lat) * sin(long); cos(lat)];
    c2 = [-sin(long); cos(long); 0];
    c3 = [-cos(lat) * cos(long); -cos(lat) * sin(long); -sin(lat)];
    R = [c1, c2, c3];

    % Rotate vector
    ref_vector_ecef = R * ref_vec_ned;

end
