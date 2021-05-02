function Q_selection(Eclipse,Q_no_eclipse,R_no_eclipse,mekf,Q_eclipse_load)

    if (Eclipse)
        % Variances
        %Q = Q_eclipse_load; % Variance of the process noise w[cycle_index]
        Q = 1e-4*diag([1 1 1 1e-3 1e-3 1e-3]);
        number_of_measurements = 6;
        % R Variances used in MEKF
        R_hat_coeff=[1e-8;1e-8;1e-8;1;1;1];
%         R_hat_coeff=10000*[.5e-3;.5e-3;.5e-3;4e-3;4e-3;4e-3;1e-3;1e-3;1e-3];
        R_hat = R_hat_coeff.*eye(number_of_measurements,number_of_measurements);
        mekf.setProcessNoiseCov(Q); %Q variance matrix
        mekf.setMeasureNoiseCov(R_hat); %R variance matrix
    else
        mekf.setProcessNoiseCov(Q_no_eclipse); %Q variance matrix
        mekf.setMeasureNoiseCov(R_no_eclipse); %R variance matrix
    end

