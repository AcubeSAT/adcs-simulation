% -----------------------------------------------------------------------------
%
%                              procedure msrCookieFinalExtended
%
%  this procedure is a struct storing all required variables
%  used for MEKF correct function.
%
%   inputs        :
%     magn_ref    - reference magnetic field vector
%     sun_ref     - reference sun position vector
%     eclipse     - existence or not of eclipse conditions
%     gyro        - gyroscope measurement
%     xsat_eci    - satellite position vector in ECI frame
%     albedo      - albedo disturbance vector
%     lamda       - Poisson distribution constant for noise
%
%   outputs       :
%     cookie      - struct including all input variables
%  ----------------------------------------------------------------------------*/

function cookie = msrCookieFinalExtended(magn_ref, sun_ref, eclipse, gyro, xsat_eci, albedo, lambda)

    cookie = struct('magn_ref', magn_ref, 'sun_ref', sun_ref, 'eclipse', eclipse, 'gyro', gyro, 'xsat_eci', xsat_eci, 'albedo', albedo, 'lambda', lambda);

end
