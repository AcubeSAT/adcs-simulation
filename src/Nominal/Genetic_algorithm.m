% =========================================================================%
%   Implementation of the Genetic Algorithm using the corresponding
%   function of MATLAB. This algorithm is used for the Gain Tuning of ADCS
%   controllers.
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
%       (fitness_function_Nominal.m)
% =========================================================================%

fit = @fitness_function_Nominal;

lb = [1 1 1 1 1 1];
ub = [100 100 100 100 100 100];
% lb = [0 0 0 0 0 0];
% ub = [1 1 1];

opts = optimoptions(@ga, ...
                    'PopulationSize', 20, ...
                    'MaxGenerations', 10, ...
                    'UseParallel', true);

[x, f] = ga(fit, 6, [], [], [], [], lb, ub, [], opts);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

% M1 = x;
% name = 'genetic.xls';
% writematrix(M1, name)