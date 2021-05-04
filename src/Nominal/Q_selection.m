function Q_selection(Eclipse,Q_no_eclipse,R_no_eclipse,mekf,Q_eclipse_load)

    if (Eclipse)
        % Variances
        %Q = Q_eclipse_load; % Variance of the process noise w[cycle_index]
        
        Q_eclipse = 8.66e-04*diag([7.68 7.68 7.68 1 1 1]);      % 6 PM
        R_eclipse = [2; 2; 2; 12; 12; 12].*eye(6,6);            % 6 PM   

%         Q_eclipse = 1e-4*diag([1 1 1 1e-3 1e-3 1e-3]);          % 11 PM
%         R_eclipse =[1e-8;1e-8;1e-8;1;1;1].*eye(6, 6);           % 11 PM
        
        mekf.setProcessNoiseCov(Q_eclipse); %Q variance matrix
        mekf.setMeasureNoiseCov(R_eclipse); %R variance matrix
    else
        mekf.setProcessNoiseCov(Q_no_eclipse); %Q variance matrix
        mekf.setMeasureNoiseCov(R_no_eclipse); %R variance matrix
    end


