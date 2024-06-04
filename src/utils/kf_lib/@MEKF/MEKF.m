%% MEKF class
% Implementation of the discrete Multiplicative Extended Kalman Filter algorithm.
% - The Jacobian of the state transition and measurement functions can either
%   be provided or approximated numerically.
% - A fading memory coefficient can be defined.
% - Enforcement of linear constraints can be imposed on the estimated parameteres.
% -


classdef MEKF < handle
    properties
        F_k % state transition function Jacobian
        H_k % measurement function Jacobian

        innov_history % Innovation history
        K % Kalman gain
        Q % process noise covariance
        R % measurement noise covariance

        a_p % fading memory coefficient
        global_state % parameters estimate
        P % parameters estimation error covariance

        % Apply projection so that:
        % A_c * global_state <= b_c
        enable_constraints
        A_c % constraint matrix
        b_c % constraints bounds

        apply_cov_sat
        global_state_sigma_min
        global_state_sigma_max
        %         eclipse

        stateTransFun_ptr % state transition function pointer
        msrFun_ptr % measurement function pointer

        stateTransFunJacob_ptr % state transition function Jacobian pointer
        msrFunJacob_ptr % measurement function Jacobian pointer

        dtheta % Parameters step size for numerical approximation of transition and measurement function Jacobian
    end

    methods

        %% Extended Kalman Filter constructor
        %  @param[in] N_params: Number of states to estimate.
        %  @param[in] N_out: Number of outputs.
        %  @param[in] stateTransFun_ptr: Pointer to state transition function. It should accept a vector (the parameters) and
        %                                  a void pointer (cookie) to pass optinally extra arguments to the state transition function.
        %  @param[in] msrFun_ptr: Pointer to measurement function. It should accept a vector (the parameters) and
        %                                  a void pointer (cookie) to pass optinally extra arguments to the measurement function.
        function this = MEKF(N_params, N_out, stateTransFun_ptr, msrFun_ptr)

            this.global_state = zeros(N_params);
            this.P = zeros(N_params, N_params);

            this.setProcessNoiseCov(eye(N_params, N_params)*1e-5);
            this.setMeasureNoiseCov(eye(N_out, N_out)*0.05);

            this.setFadingMemoryCoeff(1.0);

            this.enableParamsContraints(false);
            this.setParamsConstraints([], []);

            this.setPartDerivStep(0.001);

            this.stateTransFun_ptr = stateTransFun_ptr;
            this.msrFun_ptr = msrFun_ptr;

            this.msrFunJacob_ptr = @this.calcMsrFunJacob;
            this.stateTransFunJacob_ptr = @this.calcStateTransFunJacob;

            this.innov_history = zeros(9, 1);
        end

        %% Sets the fading memory coefficient.
        %  Note if the continuous time fading memory coefficient is 'a' then
        %  the discrete one is 'a_p = exp(a*Ts)' where Ts is the sampling period.
        %  @param[in] a_p: Fading memory coefficient.
        function setFadingMemoryCoeff(this, a_p)

            this.a_p = a_p;

        end

        %% Enables/Disables constraints in the estimation parameters.
        %  Note that the constraints must be set with 'setParamsConstraints' first.
        %  @param[in] enable_contraints: Flag that is true/false for enabling/disabling the constraints.
        function enableParamsContraints(this, enable_contraints)

            this.enable_constraints = enable_contraints;

        end

        %% Sets linear constraints in the estimation parameters.
        %  The constraints are such that D*global_state <= d
        %  @param[in] A_c: Constraint matrix.
        %  @param[in] b_c: Constraints bounds.
        function setParamsConstraints(this, A_c, b_c)

            this.A_c = A_c;
            this.b_c = b_c;

        end

        %% Sets the covariance matrix of the process noise.
        %  Note that if the continuous time process noise is R_c the discrete
        %  one is 'R = Q_c/Ts' where Ts is the sampling period.
        %  @param[in] Q: Process noise covariance matrix.
        function setProcessNoiseCov(this, Q)

            this.Q = Q;

        end

        %% Sets the covariance matrix of the measurement noise.
        %  Note that if the continuous time process noise is Q_c the discrete
        %  one is 'Q = Q_c*Ts' where Ts is the sampling period.
        %  @param[in] R: Measurement noise covariance matrix.
        function setMeasureNoiseCov(this, R)

            this.R = R;

        end

        %% Sets the state transition function Jacobian.
        %  @param[in] stateTransFunJacob_ptr: Pointer to the state transition function Jacobian. It should accept a vector (the parameters)
        %                                        and a void pointer (cookie) to pass extra arguments to the function.
        function setStateTransFunJacob(this, stateTransFunJacob_ptr)

            this.stateTransFunJacob_ptr = stateTransFunJacob_ptr;

        end

        %% Sets the measurement function Jacobian.
        %  @param[in] msrFunJacob_ptr: Pointer to the measurement function Jacobian. It should accept a vector (the parameters)
        %                                        and a void pointer (cookie) to pass extra arguments to the function.
        function setMsrFunJacob(this, msrFunJacob_ptr)

            this.msrFunJacob_ptr = msrFunJacob_ptr;

        end

        %% Sets the step for computing the partial derivatives of the state
        %  transition and/or measurement function w.r.t. the estimation parameters.
        %  @param[in] dtheta: Scalar or vector of parameters step size.
        function setPartDerivStep(this, dtheta)

            if (isscalar(dtheta))
                this.dtheta = ones(length(this.global_state), 1) * dtheta;
            else
                this.dtheta = dtheta;
            end

        end

        %% Performs the MEKF prediction (time update).
        %  @params[in] cookie: Void pointer to additional arguments needed by the
        %                      state transition function and its Jacobian (default = []).

        function predict(this, cookie, dt)

            if (nargin < 2), cookie = []; end
            % ==== Jacobian Matrix and state estimation calculation ====
            if (~isempty(cookie))
                this.F_k = this.stateTransFunJacob_ptr(this.global_state, cookie);
                this.global_state = this.stateTransFun_ptr(this.global_state, cookie);
            else
                this.F_k = this.stateTransFunJacob_ptr(this.global_state);
                this.global_state = this.stateTransFun_ptr(this.global_state);
            end
            % =====  Calculate new covariance  =====

            Phi = expm(this.F_k*dt);
            this.P = Phi * this.P * Phi' + this.Q;

        end

        %% Performs the MEKF correction (measurement update).
        %  @param[in] z: Groundtruth measurement.
        %  @params[in] cookie: Void pointer to additional arguments needed by the
        %                      measurement function and its Jacobian (default = []).

        %  This is a modified version of the classical EKF correction step.
        %  Instead of directly updating the global state (this.global_state),
        %  a local error state is calculated and then used to adjust the
        %  global state. The local error quaternion is produced using the
        %  first 3 elements of the error state.
        %  The global quaternion is updated by multiplying with the
        %  local_error_quaternion while global bias (4-6 elements of state) is updated by adding
        %  bias from the local state.

        function correct(this, z, cookie)

            if (nargin < 3), cookie = []; end

            % =====  Retrieve the measurement function Jacobian  =====

            this.H_k = this.msrFunJacob_ptr(this.global_state, cookie);
            z_hat = this.msrFun_ptr(this.global_state, cookie);

            % =====  Correction procedure =====

            Kg = this.P * this.H_k' / (this.H_k * this.P * this.H_k' + this.R); %

            z_new(1:3) = z(1:3);
            z_new(4:6) = z(7:9);

            local_error_state = Kg * (z_new' - z_hat);
            local_error_quaternion = [1; 0.5 * local_error_state(1:3)];

            this.global_state(1:4) = quatProd(this.global_state(1:4), local_error_quaternion);
            this.global_state(1:4) = this.global_state(1:4) / norm(this.global_state(1:4));
            this.global_state(5:7) = this.global_state(5:7) + local_error_state(4:6);

            N_params = 6;
            I = eye(N_params, N_params);

            this.P = (I - Kg * this.H_k) * this.P;
            this.K = Kg;

            %             % =====  Apply projection if enabled  =====
            %             proj_flag = false;
            %             D = []; % active contraints
            %             d = [];
            %             if ( this.enable_constraints && ~isempty(this.b_c) )
            %                 ind = find(this.A_c*this.global_state > this.b_c);
            %                 if (~isempty(ind))
            %                     proj_flag = true;
            %                     D = this.A_c(ind,:);
            %                     d = this.b_c(ind);
            %                 end
            %             end
            %
            %             N_params = 6;
            %             I = eye(N_params, N_params);
            %
            %             if (proj_flag)
            %                 % Kg = ( I - this.P*D'/(D*this.P*D')*D ) * Kg;
            %                 % this.global_state = this.global_state - this.P*D'/(D*this.P*D')*(D*this.global_state-d);
            %                 Kg = ( I - D'/(D*D')*D ) * Kg;
            %                 this.global_state = this.global_state - D'/(D*D')*(D*this.global_state-d);
            %             end

            %             if (proj_flag)
            %                 ind = find(this.A_c*this.global_state - this.b_c > 1e-6);
            %                 if (~isempty(ind))
            %                     warning('Some constraints were violated!');
            %                     D = this.A_c(ind,:)
            %                     d = this.b_c(ind)
            %                     global_state = this.global_state'
            %                     lambda = sort(eig(this.P))'
            %                     lambda2 = sort(svd(this.P))'
            %                     pause
            %                 end
            %             end

            %       =====  Calculate new covariance  =====
            %             this.P = (I - Kg*this.H_k) * this.P;
            %         if cookie.eclipse
            %             S_hat = zeros(9,9);
            %             window = 70;
            %             if length(this.innov_history(1,:))>=window + 1
            %                 for j=0:window-1
            %                     S_hat = S_hat + this.innov_history(:,end-j)*this.innov_history(:,end-j)';
            %                 end
            %                 S_hat = S_hat/window;
            %                 this.Q = Kg * S_hat * Kg';
            %             end
            %          end
            % end

        end

        %% Computes numerically the state transition function's Jacobian.
        %  @param[in] global_state: Parameters around which the Jacobian is computed.
        %  @params[in] cookie: Void pointer to additional arguments needed by the
        %                      state transition function's Jacobian (default = []).
        %  @return F_k: The state transition function's Jacobian.
        function F_k = calcStateTransFunJacob(this, global_state, cookie)

            if (nargin < 3), cookie = []; end

            N_params = length(global_state);
            F_k = zeros(N_params, N_params);

            dtheta_j = zeros(N_params, 1);

            for j = 1:N_params
                dtheta_j(j) = this.dtheta(j);
                Ftheta2 = this.stateTransFun_ptr(global_state+dtheta_j, cookie);
                Ftheta1 = this.stateTransFun_ptr(global_state-dtheta_j, cookie);
                F_k(:, j) = (Ftheta2 - Ftheta1) / (2 * this.dtheta(j));
                dtheta_j(j) = 0.0;
            end

        end

        %% Computes numerically the measurement function's Jacobian.
        %  @param[in] theta: Parameters around which the Jacobian is computed.
        %  @params[in] cookie: Void pointer to additional arguments needed by the
        %                      measurement function's Jacobian (default = []).
        %  @return H_k: The measurement function's Jacobian.
        function H_k = calcMsrFunJacob(this, global_state, cookie)

            if (nargin < 3), cookie = []; end

            N_params = length(global_state);
            N_out = size(this.R, 1);

            H_k = zeros(N_out, N_params);

            dtheta_j = zeros(N_params, 1);

            for j = 1:N_params
                dtheta_j(j) = this.dtheta(j);
                Htheta2 = this.msrFun_ptr(global_state+dtheta_j, cookie);
                Htheta1 = this.msrFun_ptr(global_state-dtheta_j, cookie);
                H_k(:, j) = (Htheta2 - Htheta1) / (2 * this.dtheta(j));
                dtheta_j(j) = 0.0;
            end
        end


    end
end
