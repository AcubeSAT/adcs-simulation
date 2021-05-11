function Q_selection(Eclipse,Q_no_eclipse,R_no_eclipse,mekf,Q_eclipse_load)

    if (Eclipse)
        % Variances
        %Q = Q_eclipse_load; % Variance of the process noise w[cycle_index]
        
        Q_eclipse = 1e-4*diag([1 1 1 1e-3 1e-3 1e-3]); 
        R_eclipse =[1e-7;1e-7;1e-7;1;1;1].*eye(6, 6);
        
        mekf.setProcessNoiseCov(Q_eclipse); %Q variance matrix
        mekf.setMeasureNoiseCov(R_eclipse); %R variance matrix
    else
        mekf.setProcessNoiseCov(Q_no_eclipse); %Q variance matrix
        mekf.setMeasureNoiseCov(R_no_eclipse); %R variance matrix
    end


