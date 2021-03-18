% -----------------------------------------------------------------------------
%
%                              procedure msrCookieFinal
%
%  this procedure is a struct storing all required variables for the
%  real satellite model. The satellite model is used to simulate the 
%  real attitude at any point during orbit.
%
%   inputs        :
%     magn_ref    - reference magnetic field vector
%     sun_ref     - reference sun position vector
%     eclipse     - existence or not of eclipse conditions
%     gyro        - gyroscope measurement
%
%   outputs       :
%     cookie      - struct including all input variables
%  ----------------------------------------------------------------------------*/


function cookie = msrCookieFinal(magn_ref,sun_ref,eclipse,gyro)

    cookie = struct('magn_ref',magn_ref,'sun_ref', sun_ref,'eclipse',eclipse,'gyro',gyro);
end

