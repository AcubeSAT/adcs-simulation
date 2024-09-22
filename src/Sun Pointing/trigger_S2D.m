%% ======================================================================= %%
%  Function which calculates when Detumbling is ready to be activated,
%  when in Nominal Mode.
%
%  Inputs:
%      velocity_cur          - Current Velocity of the satellite (3-axis)
%      velocity_old          - Previous Velocity of the satellite (3-axis)
%      threshold_times       - Counts the times when the threshold has been
%                              exceeded
%      threshold_exceptions  - Counts the exceptions, namely the times
%                              when the threshold is not exceeded, while
%                              in trigger mode
%      N2D_threshold         - Threshold for switching from nominal to detumbling
%       total_limit          - The total number of times the estimated angular
%                              velocity must be below the threshold before
%                              activation.
%      exceptions_limit      - The number of allowable exceptions where the
%                              estimated angular velocity exceeds the threshold
%                              but is still counted towards the total limit.
%
%  Outputs:
%      trigger_flag          - Indicates when the trigger is activated after
%                              processing
%      trigger_flag_raw      - Indicates when the trigger is activated without
%                              processing
%

%% ======================================================================= %%

function [trigger_flag, trigger_flag_raw, threshold_times, threshold_exceptions] = ...
        trigger_N2D(velocity_cur, velocity_old, threshold_times, threshold_exceptions,N2D_threshold,total_limit,exceptions_limit)

    
    

    trigger_flag = 0;
    trigger_flag_raw = 0;

    if abs(velocity_cur(1)) > N2D_threshold ...
            || abs(velocity_cur(2)) > N2D_threshold ...
            || abs(velocity_cur(3)) > N2D_threshold
        if threshold_times == 0
            threshold_times = 1;
        end
        if threshold_times >= 1
            if abs(velocity_old(1)) > N2D_threshold ...
                    || abs(velocity_old(2)) > N2D_threshold ...
                    || abs(velocity_old(3)) > N2D_threshold
                threshold_times = threshold_times + 1;
            elseif (abs(velocity_old(1)) <= N2D_threshold ...
                    && abs(velocity_old(2)) <= N2D_threshold ...
                    && abs(velocity_old(3)) <= N2D_threshold) ...
                    && threshold_exceptions < exceptions_limit
                threshold_times = threshold_times + 1;
                threshold_exceptions = threshold_exceptions + 1;
            elseif (abs(velocity_old(1)) <= N2D_threshold ...
                    && abs(velocity_old(2)) <= N2D_threshold ...
                    && abs(velocity_old(3)) <= N2D_threshold) ...
                    && threshold_exceptions >= exceptions_limit
                threshold_times = 0;
                threshold_exceptions = 0;
            end
        end
        if threshold_times >= total_limit
            trigger_flag = 1;
        end
    end

    if abs(velocity_cur(1)) > N2D_threshold ...
            || abs(velocity_cur(2)) > N2D_threshold ...
            || abs(velocity_cur(3)) > N2D_threshold

        trigger_flag_raw = 1;
    end

end







%solution 1 


% function [trigger_flag, trigger_flag_raw, threshold_times, threshold_exceptions] = ...
%         trigger_S2D(velocity_cur, velocity_old, threshold_times, threshold_exceptions,N2D_threshold,total_limit,exceptions_limit)
% 
% 
% 
%     second_threshold_times=0;
%     trigger_flag = 0;
%     trigger_flag_raw = 0;
%     margin =0.13;
% 
%       if abs(velocity_cur(1)) > N2D_threshold ...
%             || abs(velocity_cur(2)) > N2D_threshold ...
%             || abs(velocity_cur(3)) > N2D_threshold
%         if threshold_times == 0
%             threshold_times = 1;
%         end
%         if threshold_times >= 1
%             if abs(velocity_old(1)) > N2D_threshold ...
%                     || abs(velocity_old(2)) > N2D_threshold ...
%                     || abs(velocity_old(3)) > N2D_threshold
%                 threshold_times = threshold_times + 1;
%             elseif (abs(velocity_old(1)) <= N2D_threshold ...
%                     && abs(velocity_old(2)) <= N2D_threshold ...
%                     && abs(velocity_old(3)) <= N2D_threshold) ...
%                     && threshold_exceptions < exceptions_limit
%                 threshold_times = threshold_times + 1;
%                 threshold_exceptions = threshold_exceptions + 1;
%             elseif (abs(velocity_old(1)) <= N2D_threshold ... 
%                     && abs(velocity_old(2)) <= N2D_threshold ...
%                     && abs(velocity_old(3)) <= N2D_threshold) ...
%                     && threshold_exceptions >= exceptions_limit
%                 threshold_times = 0;
%                 threshold_exceptions = 0;
%             end
%         end
%       end 
%    if threshold_times >= total_limit 
% 
%        if abs(velocity_cur(1)) > margin ...
%             || abs(velocity_cur(2)) > margin ...
%             || abs(velocity_cur(3)) > margin
%         if  second_threshold_times == 0
%             second_threshold_times = 1;
%         end
%         if second_threshold_times >= 1
%             if abs(velocity_old(1)) > margin...
%                     || abs(velocity_old(2)) >margin...
%                     || abs(velocity_old(3)) > margin
%                 second_threshold_times = second_threshold_times + 1;
%             elseif (abs(velocity_old(1)) <= margin ...
%                     && abs(velocity_old(2)) <= margin ...
%                     && abs(velocity_old(3)) <= margin) ...
%                     && threshold_exceptions < exceptions_limit
%                 second_threshold_times = second_threshold_times + 1;
%                 threshold_exceptions = threshold_exceptions + 1;
%             elseif (abs(velocity_old(1)) <= margin ... 
%                     && abs(velocity_old(2)) <= margin...
%                     && abs(velocity_old(3)) <= margin) ...
%                     && threshold_exceptions >= exceptions_limit
%                 second_threshold_times=0;
%                 threshold_exceptions = 0;
%             end
%         end
% 
%        end
%        if second_threshold_times >= 10
%             trigger_flag = 1;
%        end
%    end
% 
%     if abs(velocity_cur(1)) > N2D_threshold ...
%             || abs(velocity_cur(2)) > N2D_threshold ...
%             || abs(velocity_cur(3)) > N2D_threshold
% 
%         trigger_flag_raw = 1;
%     end
% 
% 
% 
% end
% 
% 
%
%
%
% solution 2
%function [trigger_flag, trigger_flag_raw, threshold_times, threshold_exceptions] = ...
%         trigger_S2D(velocity_cur, velocity_old, threshold_times, threshold_exceptions, ...
%         N2D_threshold, total_limit, exceptions_limit)
% 
%     % Initialize variables
% 
%     trigger_flag = 0;
%     trigger_flag_raw = 0;
%     margin = 0.13;
%     persistent trend_increase_count; % Use persistent to maintain the count across function calls
% 
%     if isempty(trend_increase_count)
%         trend_increase_count = 0;
%     end
% 
%     persistent trend_decrease_count; % Use persistent to maintain the count across function calls
% 
%     if isempty(trend_decrease_count)
%         trend_decrease_count = 0;
%     end
% 
% 
% 
%       if abs(velocity_cur(1)) > margin ...
%             || abs(velocity_cur(2)) > margin ...
%             || abs(velocity_cur(3)) > margin
%           trigger_flag=1;
%       else    
%       if abs(velocity_cur(1)) > N2D_threshold ...
%             || abs(velocity_cur(2)) > N2D_threshold ...
%             || abs(velocity_cur(3)) > N2D_threshold
%         if threshold_times == 0
%             threshold_times = 1;
%         end
%         if threshold_times >= 1
%             if abs(velocity_old(1)) > N2D_threshold ...
%                     || abs(velocity_old(2)) > N2D_threshold ...
%                     || abs(velocity_old(3)) > N2D_threshold
%                 threshold_times = threshold_times + 1;
%             elseif (abs(velocity_old(1)) <= N2D_threshold ...
%                     && abs(velocity_old(2)) <= N2D_threshold ...
%                     && abs(velocity_old(3)) <= N2D_threshold) ...
%                     && threshold_exceptions < exceptions_limit
%                 threshold_times = threshold_times + 1;
%                 threshold_exceptions = threshold_exceptions + 1;
%             elseif (abs(velocity_old(1)) <= N2D_threshold ... 
%                     && abs(velocity_old(2)) <= N2D_threshold ...
%                     && abs(velocity_old(3)) <= N2D_threshold) ...
%                     && threshold_exceptions >= exceptions_limit
%                 threshold_times = 0;
%                 threshold_exceptions = 0;
%             end
%         end
%       end
% 
%     % If the threshold count exceeds the total limit, check the upper margin
%           if threshold_times >= total_limit
%           if ((abs(velocity_cur(1)) <= margin &&  abs(velocity_cur(1))>N2D_threshold))  ...
%            || ((abs(velocity_cur(2)) <= margin &&  abs(velocity_cur(2)) > N2D_threshold ))  ...
%           || (abs(velocity_cur(3)) <= margin &&  abs(velocity_cur(3))> N2D_threshold)
% 
% 
%                   % Check if velocity is increasing or decreasing
%                 if abs(velocity_cur) > abs(velocity_old)...
% 
%                     trend_increase_count = trend_increase_count + 1;
% 
% 
%                 elseif abs(velocity_cur) < abs(velocity_old)...
% 
% 
%                     trend_decrease_count=trend_decrease_count +1;
%                 end    
% 
% 
% 
% 
%            end
% 
%                    if trend_increase_count > trend_decrease_count
%                        trigger_flag=1;
%                    end
%           end
%           end
% 
%     if abs(velocity_cur(1)) > N2D_threshold || ...
%        abs(velocity_cur(2)) > N2D_threshold || ...
%        abs(velocity_cur(3)) > N2D_threshold
% 
%         trigger_flag_raw = 1;
%     end
% 
% 