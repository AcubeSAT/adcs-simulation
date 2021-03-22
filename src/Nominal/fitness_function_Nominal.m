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

function [fitness] = fitness_function_Nominal(gain_vector)
    fitness = 0;
    
    Kp_gain= 1e-05*diag([gain_vector(1) gain_vector(2) gain_vector(3)]);
    Kd_gain= 1e-04*diag([gain_vector(4) gain_vector(5) gain_vector(6)]);
 
    weight1 = 1;
    weight2 = 1;
    weight3 = 1;
    
    [instant_error_perform, Time] = Nominal_Simulation(Kp_gain, Kd_gain);
    for ut_in = 21:length(Time)
        fitness = fitness + weight1*abs(instant_error_perform(ut_in, 1)) ...
                          + weight2*abs(instant_error_perform(ut_in, 2)) ...
                          + weight3*abs(instant_error_perform(ut_in, 3));  
    end
end