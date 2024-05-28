% -----------------------------------------------------------------------------
%
%                              procedure quat_EB2OB
%
%  this procedure creates a quaternion from the ECI to Orbit frame
%
%  Inputs         :
%   q_eci_body    - Quaternion from ECI to Body frame
%   nodem         - Right Ascension of the Ascending Node
%   inclm         - Inclination
%   argpm         - Argument of Perigee
%   mm            - Mean Anomaly
%
%  Outputs        :
%   q_orbit_body  - quaternion from Orbit to Body frame
%  ----------------------------------------------------------------------------*/


function q_orbit_body = quat_EB2OB(q_eci_body, nodem, inclm, argpm, mm)

    q_orbit_eci = dcm2quat(Orbit2ECI_DCM(nodem, inclm, argpm+mm));

    q_orbit_body = quatProd(q_orbit_eci, q_eci_body);
    q_orbit_body = q_orbit_body / norm(q_orbit_body);

end
