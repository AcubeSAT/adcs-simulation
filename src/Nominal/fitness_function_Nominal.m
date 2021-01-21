function [fitness] = fitness_function_Nominal(gain_vector)
    fitness = 0;
    
    Kp_gain= 2e-05*diag([gain_vector(1) gain_vector(2) gain_vector(3)]);
    Kd_gain= 2e-03*diag([gain_vector(4) gain_vector(5) gain_vector(6)]);
    weight1 = 1;
    weight2 = 1;
    weight3 = 1;
    [instant_error_perform, Time] = Nominal_Simulation(Kp_gain, Kd_gain);
    for ut_in = 21:length(Time)
        fitness = fitness + weight1*(instant_error_perform(ut_in, 1)^2) ...
                          + weight2*(instant_error_perform(ut_in, 2)^2) ...
                          + weight3*(instant_error_perform(ut_in, 3)^2);  
    end
end