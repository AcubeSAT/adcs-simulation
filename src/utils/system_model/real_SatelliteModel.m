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
            
           
            I_inv = eye(3) / this.I;
            rw_ang_momentum_vector = [0;0;this.rw_ang_momentum];
            dydt = I_inv*(this.torq + cross(this.I*vRot + rw_ang_momentum_vector ,vRot)); 
            
            dvRot = dydt;
            Q_next = quatProd(Q,quatExp(vRot*this.dt));
            
            vRot_next = vRot + dvRot*this.dt;

            x_next = [Q_next; vRot_next];

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
            y1 = quatProd( quatInv(Q) , quatProd([0; this.magn_ref], Q) );
            
            y(1:3) = y1(2:4);
            y(1:3) = y(1:3)/norm(y(1:3));
    
            this.sun_ref = this.sun_ref/norm(this.sun_ref);

            if this.eclipse == 0
                y2 = quatProd( quatInv(Q), quatProd([0; this.sun_ref], Q) );
            
                y(7:9) = y2(2:4);
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
