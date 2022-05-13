% =========================================================================%
%   Implementation of the fitness function used for tuning R and Q matrices of Kalman filter.
%   This function is called by the Genetic Algorithm
%   (genetic_kalman.m) returning the fitness value, through which
%   the matrices are calculated.
%
%   Inputs:
%       gain_vector    - A vector where the diagonal elements of matrices
%                        R and Q are included. 
%                         
%
%   Outputs:
%       fitness         - This parameter is calculated so it can be used
%                           by the Genetic Algorithm for determining the
%                           suitable region of the desired gains. The 
%                           Genetic Algorithm "proposes" gains and the 
%                           fitness parameter indicates the fit of those 
%                           gains. The fitness parameter is the sum of
%                           the mean square errors for each euler angle, the
%                           iterations needed to converge and a parameter  
%                           that increases each time an angle deviates from stability.                             
%
% =========================================================================%

function fitness = fitness_function(x)
  
    
    % load values from sensors
    
    acc = load('accelerometer.mat');
    accelerometer = acc.accelerometer;
    gyro = load('gyroscope.mat');
    gyroscope = gyro.gyroscope;
    magn = load('magnetometer.mat');
    magnetometer = magn.magnetometer;
    
    
    % ground truth is calculated by initializing Kalman with R and Q
    % as identity matrices and calculating the mean of the last 30-40 results
    %Need to change for different accelerometer, gyro and magnetometer
    %values
    ground_x = -89.913363087529630;
    ground_y = -61.466360848669574;
    ground_z = -92.835876046292820;
    
    % Initialize Kalman 
    dt = .1;
    I = eye(3,3);
    model = SatelliteModel(dt, I);
    n_params = 7;
    number_of_measurments = 6;
    mekf = MEKF(n_params, number_of_measurments, @model.stateTransFun, @model.msrFun); % Init EKF
    mekf.global_state = [1 0 0 0 0 0 0]'; % Init state estimation
    mekf.P = eye(6,6);
    Q = diag([exp(x(7:10)) exp(x(10)) exp(x(10))]); % Init Covariance matrix
    R_hat= diag(exp(x(1:6)));
    mekf.setProcessNoiseCov(Q); % Q variance matrix
    mekf.setMeasureNoiseCov(R_hat); % R variance matrix
    mekf.setFadingMemoryCoeff(1.00); % This parameter defined the "memory" of the filter, 1 meaning regular
    mekf.setPartDerivStep(0.001);  % Arithmetical jacobian step
    % If analytical jacobian is used, define the functions
    mekf.setStateTransFunJacob(@model.stateTransFunJacob);
    mekf.setMsrFunJacob(@model.msrFunJacob);

    torq = zeros(3,1);
   
    %gyro = [1;1;1];
    accelerometer_ned = [0 0 1]';
    Eclipse = 0;
    Mag_field_ned = [1 0 0]';
    Xsat_eci = zeros(3,1);
    Albedo_inaccurate = 0;
    lambda = 0;
    state = [];
    
    % aproximate error
    error = 10;
    convergence_error = 5;
    % parameter that changes to 1 when we converge
    enter = 0;
    % returns the iterations needed to reach stability
    k = 0;
    temp = 0;
    % increases each time an angle exceeds a certain value, after stability is reached
    stability_error = 0;
    % increases each time bias exceeds a certain value 
    bias_error = 0;
    
    
    for i = 1:length(gyroscope)
        gyro = gyroscope(i,:)';
        mekf.predict(stateTransCookieFinalNominal(torq,rw_ang_momentum,gyro),dt);
        y_noise = [magnetometer(i,:);accelerometer(i,:);gyroscope(i,:)]';
        mekf.correct(y_noise, msrCookieFinalExtended(Mag_field_ned,accelerometer_ned,Eclipse,gyro,Xsat_eci,Albedo_inaccurate,lambda));
        state = [state mekf.global_state];

        temp = temp+1;
        
        euler_angles = rad2deg(quat2eul(state(1:4,:)',"XYZ"));
        euler_x = euler_angles(:,1);
        euler_y = euler_angles(:,2);
        euler_z = euler_angles(:,3);
        
        % returns the number of iterations needed to reach more stable
        % results
        if abs(euler_x(i)-ground_x)<error  &&...
           abs(euler_y(i)-ground_y)<error &&...
           abs(euler_z(i)-ground_z)<error
       
                if enter == 0
                    k = temp;
                end
                enter = 1;
                
                 
                if abs(euler_x(i)-ground_x)>convergence_error  &&...
                   abs(euler_y(i)-ground_y)>convergence_error &&...
                   abs(euler_z(i)-ground_z)>convergence_error
                    stability_error = stability_error+1;
                end
                
                if abs(state(5,i))>2  || abs(state(6,i))>2 ...
                        || abs(state(6,i))>2
                    bias_error = bias_error + 1;
                end
        end
    end
    
    % calculating mean squared error for each angle 
    mse = zeros(3);
    
    mse(1) = mean((euler_x - ground_x).^2);
    mse(2) = mean((euler_y - ground_y).^2);
    mse(3) = mean((euler_z - ground_z).^2);
    
    % calculating fitness function
    fitness = 0.3*mse(1)+3*mse(2)+0.3*mse(3) + 0.5*k + 3*stability_error +  0.5*bias_error;

end
