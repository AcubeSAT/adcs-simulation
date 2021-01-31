fit = @fitness_function_Nominal;

%lb = [55 75 25 5 5 40];
%ub = [70 95 40 15 20 55];
lb = [-1 -1 -1 -1];
ub = [1 1 1 1];

opts = optimoptions(@ga, ...
                    'PopulationSize', 20, ...
                    'MaxGenerations', 5, ...
                    'UseParallel', true);

[x, f] = ga(fit, 4, [], [], [], [], lb, ub, [], opts);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

M1 = x;
name = 'genetic.xls';
writematrix(M1, name)