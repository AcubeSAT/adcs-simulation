fit = @fitness_function_Nominal;

lb = [55 75 25 5 5 40];
ub = [70 95 40 15 20 55];

opts = optimoptions(@ga, ...
                    'PopulationSize', 20, ...
                    'MaxGenerations', 6, ...
                    'UseParallel', true);

[x, f] = ga(fit, 6, [], [], [], [], lb, ub, [], opts);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

M1 = x;
name = 'genetic.xls';
writematrix(M1, name)