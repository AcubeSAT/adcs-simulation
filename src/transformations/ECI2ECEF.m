% -----------------------------------------------------------------------------
%
%                              procedure ECI2ECEF
%
%  this procedure rotates a vector from the ECI to the ECEF frame
%
%   inputs        :
%   vec_eci       - Vector in ECI frame
%   v_vec_eci     - Velocity vector in ECI frame
%   gst           - Time in GST format
%   vel           - Velocity transfromation constant (1 = transform)
%
%   outputs       :
%   vec_ecef      - Vector in ECEF frame
%   v_vec_ecef    - Velocity vector in ECEF frame
%  ----------------------------------------------------------------------------*/


function [vec_ecef, v_vec_ecef] = ECI2ECEF(vec_eci, v_vec_eci, gst, vel)

    %% ----- Vector -----
    CGAST = cos(gst);
    SGAST = sin(gst);
    vec_ecef(1) = vec_eci(1) * CGAST + vec_eci(2) * SGAST;
    vec_ecef(2) = -vec_eci(1) * SGAST + vec_eci(2) * CGAST;
    vec_ecef(3) = vec_eci(3);

    %% ----- Velocity Vector -----
    if vel == 1
        OMEGAE = 7.29211586D-5; %Earth rotation rate in rad/s
        v_vec_ecef(1) = v_vec_eci(1) * CGAST + v_vec_eci(2) * SGAST + vel * OMEGAE * vec_ecef(2);
        v_vec_ecef(2) = -v_vec_eci(1) * SGAST + v_vec_eci(2) * CGAST - vel * OMEGAE * vec_ecef(1);
        v_vec_ecef(3) = v_vec_eci(3);
    else
        v_vec_ecef = zeros(3, 1);
    end
    
end
