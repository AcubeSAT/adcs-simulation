%% Satellite model
%

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

            % if (nargin < 3), cookie = []; end
            
            % Local Torque parameter
            this.torq = cookie.torq;
            this.rw_ang_momentum=cookie.rw_ang_momentum;
            % Local quaternion and angular velocity parameters
            Q = x(1:4);
            vRot = x(5:7);
            
            
            % Calculating angular velocity dot based on dynamics model
%             dvRot = this.I\(this.torq + cross(this.I*vRot,vRot));
%             dydt(5) = (1 / this.I(1, 1)) * (-(this.I(3, 3) - this.I(2, 2)) * vRot(2) * vRot(3) + this.torq(1));
%             dydt(6) = (1 / this.I(2, 2)) * (-(this.I(1, 1) - this.I(3, 3)) * vRot(1) * vRot(3) + this.torq(2));
%             dydt(7) = (1 / this.I(3, 3)) * (-(this.I(2, 2) - this.I(1, 1)) * vRot(1) * vRot(2) + this.torq(3));
            I_inv = eye(3) / this.I;
            rw_ang_momentum_vector = [0;0;this.rw_ang_momentum];
            dydt = I_inv*(this.torq + cross(this.I*vRot + rw_ang_momentum_vector ,vRot)); 
            
            dvRot = dydt;
%             dvRot
            % Calculating quaternion dot based on kinematics model
            Q_next = quatProd(Q,quatExp(vRot*this.dt));
            
            % Calculating next quaterniong and angular velocity
            vRot_next = vRot + dvRot*this.dt;

            x_next = [Q_next; vRot_next];

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
            
            Q = x(1:4); %quaternion
            vRot = x(5:7); %angular velocity

            y = zeros(6,1); %measurements
            y(4:6) = vRot; %gyro measurements should be equal to angular velocity [wrong]
            
            % Calculating magnetometer measurements
            this.magn_ref=this.magn_ref/norm(this.magn_ref);
            y1 = quatProd( quatInv(Q) , quatProd([0; this.magn_ref], Q) );
            
            y(1:3) = y1(2:4); % quat2dcm(Q')*this.magn_ref;
            y(1:3) = y(1:3)/norm(y(1:3));
            
            %This code part models the sun sensor FOV, but it is not
            %finished yet. At the moment in determination we donot nadir
            %point hence sun sensor takes little values. This
            %implementation should be revisited when Determination is
            %integrated with control
            this.sun_ref=this.sun_ref/norm(this.sun_ref);
            Q_sun_sensor = [Q(1),Q(2),-Q(3),-Q(4)];
            quat_fov=quatProd(quatconj(Q_sun_sensor),[0; this.sun_ref]);      
            fov_angle = 2*asin(norm(quat_fov(2:4)))*180/pi;
            
            if(abs(fov_angle) < 60)
                 inside_fov = true;
            else
                inside_fov = false;
            end

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
        sun_ref
        eclipse
        rw_ang_momentum
    end
    
end
