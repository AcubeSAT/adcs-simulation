fit = @fitness_function_Nominal;

lb = [45 55 15 0 0 35];
ub = [85 100 50 25 30 90];

opts = optimoptions(@ga, ...
                    'PopulationSize', 20, ...
                    'MaxGenerations', 5, ...
                    'UseParallel', true);

[x, f] = ga(fit, 6, [], [], [], [], lb, ub, [], opts);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

M1 = x;
name = 'genetic.xls';
writematrix(M1, name)