fit = @fitness_function_Nominal;

lb = [0 0 0 0 0 0];
ub = [100 100 100 100 100 100];

opts = optimoptions(@ga, ...
                    'PopulationSize', 20, ...
                    'MaxGenerations', 4, ...
                    'UseParallel', true);

[x, f] = ga(fit, 6, [], [], [], [], lb, ub, [], opts);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

M1 = x;
name = 'genetic.xls';
writematrix(M1, name)