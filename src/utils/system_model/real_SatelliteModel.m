% ======================================================================== %
%     Class that implements the real groundtruth space environment
% 
%     Properties:
%         dt                 - Simulation time step
%         I                  - Inertia
%         magn_ref           - Magnetic field
%         torq               - Applied torque
%         sun_ref            - Sun position
%         eclipse            - Existence or not of eclipse conditions
%         rw_ang_momentum    - RW angular momentum
% ======================================================================== %

classdef real_SatelliteModel

    methods (Access = public)
        
        %% Constructor.
        % @param[in] dt: sampling time.
        function this = real_SatelliteModel(dt, I)
            
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
            this.rw_ang_momentum = cookie.rw_ang_momentum;
            Q = x(1:4);
            vRot = x(5:7);
            
            rw_ang_momentum_vector = [0;0;this.rw_ang_momentum];

            
            % -------- Runge Kutta 4th Order  -----------------------------
            % First step
            vRot_quat = [0 vRot(1) vRot(2) vRot(3)];
            k1_v = this.I \ (this.torq - cross(vRot, (this.I * vRot + rw_ang_momentum_vector)));
            k1_q = 0.5 * quatProd(Q, vRot_quat);

            % Second step
            vRot2 = vRot + 0.5 * this.dt * k1_v;
            vRot2_quat = [0 vRot2(1) vRot2(2) vRot2(3)];
            q2 = Q + 0.5 * this.dt * k1_q;
            k2_v = this.I \ (this.torq - cross(vRot2, (this.I * vRot2 + rw_ang_momentum_vector)));
            k2_q = 0.5 * quatProd(q2, vRot2_quat);

            % Third step
            vRot3 = vRot + 0.5 * this.dt * k2_v;
            vRot3_quat = [0 vRot3(1) vRot3(2) vRot3(3)];
            q3 = Q + 0.5 * this.dt * k2_q;
            k3_v = this.I \ (this.torq - cross(vRot3, (this.I * vRot3 + rw_ang_momentum_vector)));
            k3_q = 0.5 * quatProd(q3, vRot3_quat);

            % Fourth step
            vRot4 = vRot + this.dt * k3_v;
            vRot4_quat = [0 vRot4(1) vRot4(2) vRot4(3)];
            q4 = Q + this.dt * k3_q;
            k4_v = this.I \ (this.torq - cross(vRot4, (this.I * vRot4 + rw_ang_momentum_vector)));
            k4_q = 0.5 * quatProd(q4, vRot4_quat);

            % Combine steps
            new_vRot = vRot + (this.dt / 6.0) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
            new_q = Q + (this.dt / 6.0) * (k1_q + 2 * k2_q + 2 * k3_q + k4_q);
            new_q = new_q / norm(new_q);

            x_next = [new_q; new_vRot];
            %--------------------------------------------------------------

          

            % ------- Euler integration method ----------------------------
            % I_inv = eye(3) / this.I;
            %
            % dydt = I_inv*(this.torq + cross(this.I*vRot + rw_ang_momentum_vector ,vRot)); 
            % 
            % dvRot = dydt;
            % vRot_next = vRot + dvRot*this.dt;
            % 
            % Q_next = quatProd(Q,quatExp(vRot*this.dt));      
            % x_next = [Q_next; vRot_next];
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
            
            Q = x(1:4); 
            vRot = x(5:7); 

            y = zeros(6,1); 
            y(4:6) = vRot;
            
            this.magn_ref=this.magn_ref/norm(this.magn_ref);
            
            y(1:3) = rotate_vector(Q,this.magn_ref);
            
            y(1:3) = y(1:3)/norm(y(1:3));
    
            this.sun_ref = this.sun_ref/norm(this.sun_ref);

            if this.eclipse == 0
                y(7:9) = rotate_vector(Q,this.sun_ref);
       
                y(7:9) = y(7:9)/norm(y(7:9));
            else
                y(7:9) = zeros(3,1);
            end

        end


        %% Measurement function Jacobian.
        % @param[in] x: state at time k
        % @param[in] cookie: optional arguments.
        % @param[out] H: measurement function Jacobian.
        %
        function H = msrFunJacob(this, x, cookie)
            
            error('[SatelliteModel::msrFunJacob]: Not implemented!');
            
        end


        %% State transition function Jacobian.
        % @param[in] x: state at time k
        % @param[in] cookie: optional arguments.
        % @param[out] F: state transition function Jacobian.
        %
        function F = stateTransFunJacob(this, x, cookie)
            
            error('[SatelliteModel::msrFunJacob]: Not implemented!');
            
        end

 
    end

    properties (Access = protected)
        
        dt
        I
        magn_ref
        torq
        sun_ref
        eclipse
        rw_ang_momentum
    end
    
end
