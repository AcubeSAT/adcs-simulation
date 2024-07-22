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