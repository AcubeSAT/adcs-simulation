% =========================================================================%
%   Implementation of the Genetic Algorithm using the corresponding
%   function of MATLAB. This algorithm is used for the Gain Tuning of
%   matrices Q and R of the Kalman Filter.
%
%   Inputs:
%       The options for the each run, namely:
%           Population size     - Indicates the size of each generation's
%                                    population
%           Max Generations     - Indicates the maximum number of
%                                   generations
%           Use Parallel        - Set to TRUE if we want parellel
%                                   computations
%   Ouputs:
%       The gain vector which includes the desired gains
%       The result of the fitness parameter in fitness function
%       (fitness_function.m)
% =========================================================================%


fit = @fitness_function;


%Lower and upper limits
lb = ones(10, 1) * (-10);
ub = ones(10, 1) * 10;

opts = optimoptions(@ga, ...
    'PopulationSize', 12, ...
    'MaxGenerations', 2, ...
    'UseParallel', true);

[gain_vector, f] = ga(fit, 10, [], [], [], [], lb, ub, [], opts);
