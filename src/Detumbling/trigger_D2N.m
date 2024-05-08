%% ======================================================================= %%
 %  Function which calculates when Nadir Pointing is ready to be activated,
 %  when in Detumbling Mode.
 %  
 %  Inputs: 
 %      threshold_times       - Counts the times when the threshold has been
 %                              exceeded
 %      threshold_exceptions  - Counts the exceptions, namely the times
 %                              when the threshold is not exceeded, while 
 %      w_b_ob                - Angular rate of satellite relative to ECI frame
 %      w_b_ob_Bdot           - Angular velocity estimation using Bdot  
 %      w_b_ob_Bdot_previous  - Angular velocity estimation using Bdot at
 %                              previous timestep
 %      D2N_threshold         - Threshold angular velocity [rad/sec]                    
 %
 %  Outputs:
 %      bdot_activation       - Indicates when the trigger is activated,
 %                              based on the estimated angular velocity
 %      nonBdot_activation    - Indicates when the trigger is
 %                              activated,based on the real angular velocity
 %
 %% ======================================================================= %%

function [bdot_activation,nonBdot_activation,threshold_times,threshold_exceptions] = trigger_D2N(threshold_times,threshold_exceptions,w_b_ob,w_b_ob_Bdot,w_b_ob_Bdot_previous,D2N_threshold)
        bdot_activation = 0;
        nonBdot_activation = 0;
        total_limit = 520;
        exceptions_limit = 30;
        
            if abs(w_b_ob_Bdot(1)) < D2N_threshold ... 
                    && abs(w_b_ob_Bdot(2)) < D2N_threshold ... 
                        && abs(w_b_ob_Bdot(3)) < D2N_threshold 
                 if threshold_times == 0 
                    threshold_times = 1; 
                 end 
                 if threshold_times >= 1    
                     if abs(w_b_ob_Bdot_previous(1)) < D2N_threshold ... 
                        && abs(w_b_ob_Bdot_previous(2)) < D2N_threshold ... 
                            && abs(w_b_ob_Bdot_previous(3)) < D2N_threshold 
                        threshold_times = threshold_times + 1;
                     elseif (abs(w_b_ob_Bdot_previous(1)) >= D2N_threshold ... 
                             || abs(w_b_ob_Bdot_previous(2)) >= D2N_threshold ... 
                                || abs(w_b_ob_Bdot_previous(3)) >= D2N_threshold) ...
                                    && threshold_exceptions < exceptions_limit
                        threshold_times = threshold_times + 1;
                        threshold_exceptions = threshold_exceptions + 1;
                     elseif (abs(w_b_ob_Bdot_previous(1)) >= D2N_threshold ... 
                             || abs(w_b_ob_Bdot_previous(2)) >= D2N_threshold ... 
                                || abs(w_b_ob_Bdot_previous(3)) >= D2N_threshold) ...
                                    && threshold_exceptions >= exceptions_limit
                        threshold_times = 0;
                        threshold_exceptions = 0; 
                     end 
                 end 
                 if threshold_times >= total_limit
                     bdot_activation = 1;
                 end     
            end 
            if abs(w_b_ob(1)) < D2N_threshold ... 
                    && abs(w_b_ob(2)) < D2N_threshold ... 
                        && abs(w_b_ob(3)) < D2N_threshold 
                    
                    nonBdot_activation = 1;
            end

end

