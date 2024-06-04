% =========================================================================%
%   Implementation of the fitness function used for the Nominal Simulation.
%   This function is called by the Genetic Algorithm
%   (Genetic_algorithm.m) returning the fitness parameter, through which
%   the controller's gains are calculated.
%
%   Inputs:
%       gain_vector     - A vector where the gains until final calculation
%                        are included and changed by the Genetic Algorithm.
%
%   Outputs:
%       fitness         - This parameter is calculated so it can be used
%                           by the Genetic Algorithm for determining the
%                           suitable region of the desired gains. The
%                           Genetic Algorithm "proposes" gains and the
%                           fitness parameter indicates the fit of those
%                           gains. This parameter is formed based on our
%                           definition on "desired gains"
%                           (ex. the fitness parameter is the sum of
%                           all Absolute Performance Instant Errors. The
%                           larger the value of the fitness function, the
%                           worst the region of the gains)
%
% =========================================================================%

function [fitness] = fitness_function_Nadir_Pointing(gain_vector)
    fitness = 0;

    % Kp_gain= 1e-05*diag([gain_vector(1) gain_vector(2) gain_vector(3)]);
    % Kd_gain= 1e-04*diag([gain_vector(4) gain_vector(5) gain_vector(6)]);

    Kp_gain = 1e-04 * diag([3, 5, 6]);
    Kd_gain = 1e-04 * diag([9.6, 47.6, 70]);

    n_dim_error = 6;
    number_of_measurments = 6;

    Q_no_eclipse = gain_vector(1) * eye(n_dim_error, n_dim_error);
    R_no_eclipse = [gain_vector(2); gain_vector(2); gain_vector(2); gain_vector(3); gain_vector(3); gain_vector(3)] .* eye(number_of_measurments, number_of_measurments);
    Q_eclipse_load = gain_vector(4) * eye(n_dim_error, n_dim_error);
    R_hat = [gain_vector(5); gain_vector(5); gain_vector(5); gain_vector(6); gain_vector(6); gain_vector(6)] .* eye(number_of_measurments, number_of_measurments);

    [instant_error_perform, Time] = Nadir_Pointing_function(Kp_gain, Kd_gain, Q_no_eclipse, R_no_eclipse, Q_eclipse_load, R_hat);
    for ut_in = 21:length(Time)

        if abs(instant_error_perform(ut_in, 1)) >= 100 ...
                || abs(instant_error_perform(ut_in, 2)) >= 100 ...
                || abs(instant_error_perform(ut_in, 3)) >= 100
            weight1 = 120;
            weight2 = 120;
            weight3 = 120;

        elseif abs(instant_error_perform(ut_in, 1)) >= 70 ...
                || abs(instant_error_perform(ut_in, 2)) >= 70 ...
                || abs(instant_error_perform(ut_in, 3)) >= 70

            weight1 = 24;
            weight2 = 24;
            weight3 = 24;

        elseif abs(instant_error_perform(ut_in, 1)) >= 50 ...
                || abs(instant_error_perform(ut_in, 2)) >= 50 ...
                || abs(instant_error_perform(ut_in, 3)) >= 50

            weight1 = 6;
            weight2 = 6;
            weight3 = 6;

        elseif abs(instant_error_perform(ut_in, 1)) >= 30 ...
                || abs(instant_error_perform(ut_in, 2)) >= 30 ...
                || abs(instant_error_perform(ut_in, 3)) >= 30

            weight1 = 2;
            weight2 = 2;
            weight3 = 2;

        elseif abs(instant_error_perform(ut_in, 1)) >= 10 ...
                || abs(instant_error_perform(ut_in, 2)) >= 10 ...
                || abs(instant_error_perform(ut_in, 3)) >= 10

            weight1 = 1;
            weight2 = 1;
            weight3 = 1;

        else
            weight1 = 0.1;
            weight2 = 0.1;
            weight3 = 0.1;
        end

        fitness = fitness + weight1 * abs(instant_error_perform(ut_in, 1)) ...
                          + weight2 * abs(instant_error_perform(ut_in, 2)) ...
                          + weight3 * abs(instant_error_perform(ut_in, 3));


    end
end