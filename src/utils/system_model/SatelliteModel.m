% ======================================================================== %
%     Class that implements the space environment used by the satellite
%
%     Properties:
%         dt                 - Simulation time step
%         I                  - Inertia
%         magn_ref           - Magnetic field
%         torq               - Applied torque
%         gyro               - Gyro measurements
%         sun_ref            - Sun position
%         eclipse            - Existence or not of eclipse conditions
% ======================================================================== %

classdef SatelliteModel

    methods (Access = public)

        %% Constructor.
        % @param[in] dt: sampling time.
        function this = SatelliteModel(dt, I)

            this.dt = dt;
            this.I = I;

        end

        %% State transition function.
        % @param[in] x: state at time k.
        % @param[in] cookie: optional arguments.
        % @param[out] x_next: state at time k+1.
        %
        function x_next = stateTransFun(this, x, cookie)

            this.torq = cookie.torq;
            this.gyro = cookie.gyro;

            bias = x(5:7);
            Q = x(1:4);
            vRot = this.gyro - bias;
            
            
            % -------- Runge Kutta 4th Order  -----------------------------
            vRot_quat = [0 vRot(1) vRot(2) vRot(3)];

            k1_q = 0.5 * quatProd(Q,vRot_quat);
            k2_q = 0.5 * quatProd(Q + 0.5 * this.dt * k1_q, vRot_quat);
            k3_q = 0.5 * quatProd(Q + 0.5 * this.dt * k2_q, vRot_quat);
            k4_q = 0.5 * quatProd(Q + this.dt * k3_q, vRot_quat);

            q_new = Q + (this.dt / 6.0) * (k1_q + 2 * k2_q + 2 * k3_q + k4_q);
            Q_next = q_new / norm(q_new);

            x_next = [Q_next; bias];
            % -------------------------------------------------------------

           

            % ------- Euler integration method ----------------------------
            % Q_next = quatProd( Q, quatExp(vRot*this.dt));
            % x_next = [Q_next; bias];
            % -------------------------------------------------------------

        end

        %% Measurement function.
        % @param[in] x: state at time k.
        % @param[in] cookie: optional arguments.
        % @param[out] y: measurements at time k.
        %
        function y = msrFun(this, x, cookie)

            this.magn_ref = cookie.magn_ref;
            this.sun_ref = cookie.sun_ref;
            this.eclipse = cookie.eclipse;
            this.magn_ref = this.magn_ref / norm(this.magn_ref);
            this.sun_ref = this.sun_ref / norm(this.sun_ref);
            this.gyro = cookie.gyro;

            Q = x(1:4);

            y = zeros(6, 1);

            y(1:3) = rotate_vector(Q, this.magn_ref);
            y(1:3) = y(1:3) / norm(y(1:3));

            if this.eclipse == 0
                y(4:6) = css_compensation(this.sun_ref, Q, cookie.xsat_eci, cookie.albedo, cookie.lambda);
            else
                y(4:6) = zeros(3, 1);
            end

        end

        %% Measurement function Jacobian.
        % @param[in] x: state at time k
        % @param[in] cookie: optional arguments.
        % @param[out] H: measurement function Jacobian.
        %
        function H = msrFunJacob(this, x, cookie)

            z_hat = this.msrFun(x, cookie);

            H = [skew(z_hat(1:3)), zeros(3, 3); ...
                skew(z_hat(4:6)), zeros(3, 3)];
        end

        %% State transition function Jacobian.
        % @param[in] x: state at time k
        % @param[in] cookie: optional arguments.
        % @param[out] F: state transition function Jacobian.
        %
        function F = stateTransFunJacob(this, x, cookie)

            F = [-skew(cookie.gyro-x(5:7)), -eye(3, 3); zeros(3, 6)];

        end


    end

    properties (Access = protected)

        dt
        I
        magn_ref
        torq
        gyro
        sun_ref
        eclipse
    end

end
