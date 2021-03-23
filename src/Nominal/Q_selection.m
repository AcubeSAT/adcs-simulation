function [exited_eclipse,entered_eclipse] = Q_selection(Eclipse,Q_no_eclipse,R_no_eclipse,mekf,Q_eclipse_load,exited_eclipse,entered_eclipse,set_Q_selection_method)


if set_Q_selection_method == "IdealQ"   
    if (Eclipse)
        % Variances
        Q = Q_eclipse_load; % Variance of the process noise w[cycle_index]
        number_of_measurments = 9;
        % R Variances used in MEKF
        % R_hat_coeff=[1e-3;1e-3;1e-3;8e-3;8e-3;8e-3;5e-3;5e-3;5e-3];
        R_hat_coeff=10000*[.5e-3;.5e-3;.5e-3;4e-3;4e-3;4e-3;1e-3;1e-3;1e-3];
        R_hat = R_hat_coeff.*eye(number_of_measurments,number_of_measurments);
        mekf.setProcessNoiseCov(Q); %Q variance matrix
        mekf.setMeasureNoiseCov(R_hat); %R variance matrix
    else
        mekf.setProcessNoiseCov(Q_no_eclipse); %Q variance matrix
        mekf.setMeasureNoiseCov(R_no_eclipse); %R variance matrix
    end
elseif set_Q_selection_method == "AdaptiveQ"
    number_of_measurments = 9;
    if (entered_eclipse == false) && (Eclipse) ~= 0             
        Q_eclipse =  eye(6,6);
        Q_eclipse(4:6,4:6) = 0.5e-06*eye(3,3);

        R_hat_coeff=10000*[.5e-3;.5e-3;.5e-3;4e-3;4e-3;4e-3;1e-3;1e-3;1e-3];
        R_hat = R_hat_coeff.*eye(number_of_measurments,number_of_measurments);
                  
        mekf.setProcessNoiseCov(Q_eclipse); %Q variance matrix
        mekf.setMeasureNoiseCov(R_hat); %R variance matrix
        entered_eclipse = true;
    end
    if (exited_eclipse == false) && (entered_eclipse == true) && (Eclipse) == 0
        Q_outside =  0.5e-05*eye(6,6);
        Q_outside(4:6,4:6) = 0.5e-07*eye(3,3);

        mekf.setProcessNoiseCov(Q_outside);
        mekf.setMeasureNoiseCov(R_no_eclipse); %R variance matrix
        exited_eclipse = true;
        entered_eclipse = false;
    end
end

