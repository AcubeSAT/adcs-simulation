fit = @fitness_function_Nominal;

lb = [0 0 0 0 0 0];
ub = [100 100 100 100 100 100];

opts = optimoptions(@ga, ...
                    'PopulationSize', 12, ...
                    'MaxGenerations', 12, ...
                    'UseParallel', true);

[x, f] = ga(fit, 6, [], [], [], [], lb, ub, [], opts);

