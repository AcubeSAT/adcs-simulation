function [fitness] = fitness_function_Nominal(Q_vector)
    fitness = 0;
    
    %Kp_gain= 2e-05*diag([gain_vector(1) gain_vector(2) gain_vector(3)]);
    %Kd_gain= 2e-03*diag([gain_vector(4) gain_vector(5) gain_vector(6)]);
    Q0 = [Q_vector(1) Q_vector(2) Q_vector(3) Q_vector(4)];
    weight1 = 1;
    weight2 = 1;
    weight3 = 1;
    Kp_gain =2e-05*diag([62.75,82.95,32]);%!!!% 6 -> good
    Kd_gain=2e-03*diag([10.3,14.7,48.65]);%!!!%
    %[instant_error_perform, Time] = Nominal_Simulation(Kp_gain, Kd_gain);
    [instant_error_perform, Time] = Nominal_Simulation(Kp_gain, Kd_gain, Q0');
    for ut_in = 21:length(Time)
        fitness = fitness + weight1*abs(instant_error_perform(ut_in, 1)) ...
                          + weight2*abs(instant_error_perform(ut_in, 2)) ...
                          + weight3*abs(instant_error_perform(ut_in, 3));  
    end
end