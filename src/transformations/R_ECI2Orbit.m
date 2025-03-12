% -----------------------------------------------------------------------------
%
%                              procedure ECI2Orbit
%
%  this procedure rotates a vector from the ECI to the Orbit frame
%
%  Inputs            :
%   ref_vector_eci   - Vector in ECI frame
%   nodeo            - Right Ascension of the Ascending Node
%   incl             - Inclination
%   argp_mo          - Addition of Argument of Perigee with Mean Anomaly
%
%  Outputs           :
%   ref_vector_orbit - Vector in Orbit frame
%
%
%  Based on AOCS_DDJF
%  ----------------------------------------------------------------------------*/

function R_IO = R_ECI2Orbit(nodeo, incl, argp_mo)

    % Construct rotation matrix
    R = zeros(3, 3);
    R(1, 1) = -sin(argp_mo) * sin(nodeo) * cos(incl) + cos(argp_mo) * cos(nodeo);
    R(1, 2) = sin(argp_mo) * cos(incl) * cos(nodeo) + sin(nodeo) * cos(argp_mo);
    R(1, 3) = sin(incl) * sin(argp_mo);

    R(2, 1) = -sin(argp_mo) * cos(nodeo) - sin(nodeo) * cos(incl) * cos(argp_mo);
    R(2, 2) = -sin(argp_mo) * sin(nodeo) + cos(incl) * cos(argp_mo) * cos(nodeo);
    R(2, 3) = sin(incl) * cos(argp_mo);

    R(3, 1) = sin(incl) * sin(nodeo);
    R(3, 2) = -sin(incl) * cos(nodeo);
    R(3, 3) = cos(incl);

    % Adjust rotation matrix
    R_IO = [-R(1, :); R(3, :); R(2, :)];

end
