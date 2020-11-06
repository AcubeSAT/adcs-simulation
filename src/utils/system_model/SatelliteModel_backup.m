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
            
             % Bias propagation bdot = 0
            bias = x(8:10);
            
            % Local quaternion and angular velocity parameters
            Q = x(1:4);
            vRot = x(5:7) - bias; % correction with bias
                        
            % Local quaternion dot and angular velocity dot parameter
            dydt = zeros(7,1);
            
% dydt(2) = (1 / (2 * norm(Q))) * (vRot(1) * Q(1) - vRot(2) * Q(4) + vRot(3) * Q(3));
% dydt(3) = (1 / (2 * norm(Q))) * (vRot(1) * Q(4) + vRot(2) * Q(1) - vRot(3) * Q(2));
% dydt(4) = (1 / (2 * norm(Q))) * (-vRot(1) * Q(3) + vRot(2) * Q(2) + vRot(3) * Q(1));
% dydt(1) = (1 / (2 * norm(Q))) * (-vRot(1) * Q(2) - vRot(2) * Q(3) - vRot(3) * Q(4));
% 
% Q_next=Q+dydt(1:4)*this.dt;

%             dQ = 0.5*quatProd([0; vRot], Q);
%             Q_next=Q+dQ*this.dt;
            % Calculating angular velocity dot based on dynamics model
           dydt(5) = (1 / this.I(1, 1)) * (-(this.I(3, 3) - this.I(2, 2)) * vRot(2) * vRot(3) + this.torq(1));
           dydt(6) = (1 / this.I(2, 2)) * (-(this.I(1, 1) - this.I(3, 3)) * vRot(1) * vRot(3) + this.torq(2));
           dydt(7) = (1 / this.I(3, 3)) * (-(this.I(2, 2) - this.I(1, 1)) * vRot(1) * vRot(2) + this.torq(3));
            
            
            dvRot = dydt(5:7);

%             vRot2 = 2*quatProd(dQ, quatInv(Q));
%             vRot2 = vRot2(2:4);

%             Q_next = quatProd( quatExp(vRot2*this.dt), Q );
            % Calculating quaternion dot based on kinematics model
            Q_next = quatProd( quatExp(vRot*this.dt), Q );
            
            % Calculating next quaterniong and angular velocity
            vRot_next = vRot + dvRot*this.dt;

            x_next = [Q_next; vRot_next; bias];

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
            y1 = quatProd( quatInv(Q), quatProd([0; this.magn_ref], Q) );
            
            y(1:3) = y1(2:4); % quat2dcm(Q')*this.magn_ref;
            y(1:3) = y(1:3)/norm(y(1:3));
            
            %This code part models the sun sensor FOV, but it is not
            %finished yet. At the moment in determination we donot nadir
            %point hence sun sensor takes little values. This
            %implementation should be revisited when Determination is
            %integrated with control
            quat_fov=quatProd(quatInv(Q),[0; this.sun_ref]);
            quat_fov=quat_fov/norm(quat_fov);
            fov_angle = 2*atan2(norm(quat_fov(2:4)),quat_fov(1))*180/pi;
%             abs(fov_angle-180) < 60

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
        
    end
    
end
