function [bdot_activation,nonBdot_activation,threshold_times,threshold_exceptions] = D2N_trigger(threshold_times,threshold_exceptions,w_b_ob,w_b_ob_Bdot,w_b_ob_Bdot_previous,D2N_threshold)
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

