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