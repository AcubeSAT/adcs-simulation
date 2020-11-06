%% Satellite model
%

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

            % if (nargin < 3), cookie = []; end
            
            % Local Torque parameter
            this.torq = cookie.torq;
            this.gyro = cookie.gyro;
            
            % Bias propagation bdot = 0
            bias = x(5:7);
            
            % Local quaternion and angular velocity parameters
            Q = x(1:4);
            vRot = this.gyro - bias; % correction with bias
                                
            % Calculating quaternion dot based on kinematics model
            Q_next = quatProd( Q, quatExp(vRot*this.dt));
            
            x_next = [Q_next; bias];

        end


        %% Measurement function.
        % @param[in] x: state at time k.
        % @param[in] cookie: optional arguments.
        % @param[out] y: measurements at time k.
        %
        function y = msrFun(this, x, cookie)

            % if (nargin < 3), cookie = []; end
            
            % Initializing all needed parameters
            this.magn_ref = cookie.magn_ref; %magnetic field
            this.sun_ref = cookie.sun_ref; %sun position
            this.eclipse = cookie.eclipse; %eclipse calulation
            this.magn_ref=this.magn_ref/norm(this.magn_ref);
            this.sun_ref=this.sun_ref/norm(this.sun_ref);
            this.gyro = cookie.gyro;
            
            Q = x(1:4); %quaternion
            
            y = zeros(6,1); %measurements
             
            %y(4:6) = this.gyro - x(5:7); %      Should I update?   We
            %don't propagate angular velocity
            y(4:6) = 0;   
            
            % Calculating magnetometer measurements
            y1 = quatProd( quatInv(Q), quatProd([0; this.magn_ref], Q) );
            
            y(1:3) = y1(2:4); % quat2dcm(Q')*this.magn_ref;
            y(1:3) = y(1:3)/norm(y(1:3));
            

            %If we have eclipse, we dont have a sun estimate
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
        
        dt % simulation time step
        I % inertia
        magn_ref
        torq
        gyro
        sun_ref
        eclipse
        rw_ang_momentum
    end
    
end
