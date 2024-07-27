% -----------------------------------------------------------------------------
%
%                              procedure stateTransCookieFinalNominal
%
%  this procedure is a struct storing all required inputs
%  for the satellite state model. The satellite model is used
%  for the MEKF predict function and to simulate the real attitude
%  at any point during orbit.
%
%   inputs        :
%     torq        - commanded torque vector
%     rw_ang_momentum - reaction wheel angular momentum
%     gyro        - gyroscope measurement
%
%   outputs       :
%     cookie      - struct including all input variables
%  ----------------------------------------------------------------------------*/

function cookie = stateTransCookieFinalNominal(torq, rw_ang_momentum, gyro)

    cookie = struct('torq', torq, 'rw_ang_momentum', rw_ang_momentum, 'gyro', gyro);

end
