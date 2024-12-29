
function [Bdot_body,torq,V_mtq, I_mtq, P_thermal_mtq,M] = ...
        Control_Thomson(w_b_ib, B2,B, mtq_max)

   
Kd=100;
Ks=1;
w_ref=0.087;
 Bx=acos(B(1)/norm(B));
 Bx2=acos(B2(1)/norm(B2));

    By=acos(B(2)/norm(B));
    By2=acos(B2(2)/norm(B2));
    
    Bz=acos(B(3)/norm(B));
    Bz2=acos(B2(3)/norm(B2));

    Bdot_body = ([Bx2;By2;Bz2] - [Bx;By;Bz]) / 0.1;
  

M(3,1)= Kd * Bdot_body(3);
         if (abs(B(2))>abs(B(1)))
            M(1,1)= Ks*(abs(w_b_ib(3))-w_ref)*sign(B(2));
            M(2,1)= 0;
         elseif (abs(B(2))<abs(B(1)))
            M(1,1)= 0;
            M(2,1)=-Ks*(abs(w_b_ib(3))-w_ref)*sign(B(1));
         end  

         M = mtq_scaling(M, mtq_max);
         T_magnetic_effective = cross(M, B);
         torq = T_magnetic_effective;

             %% Calculate V, I, P of MTQ's
    [V_mtq, I_mtq, P_thermal_mtq] = mtq_model(M);

  %% Angle between Z-Axis of b.f. and Z-Axis of o.f.
        % x_ob=R_OB(:,1); 
        % y_ob=R_OB(:,2);
        % z_ob=R_OB(:,3); % Extract the third column of R_OB(z-axis of orbit in body frame)
        % z_body=[0; 0; 1]; % z-axis of the body frame
        % cos_theta_x=dot(x_ob,z_body)/(norm(x_ob)*norm(z_body)); % Calculate cosine of the angle
        % cos_theta_y=dot(y_ob,z_body)/(norm(y_ob)*norm(z_body));
        % cos_theta_z=dot(z_ob,z_body)/(norm(z_ob)*norm(z_body));
        % 
        % theta_x=acos(cos_theta_x); % compute angle in radians
        % theta_y=acos(cos_theta_y);
        % theta_z=acos(cos_theta_z);
        % theta_deg_x=rad2deg(theta_x); %Convert to degrees
        % theta_deg_y=rad2deg(theta_y) ;
        % theta_deg_z=rad2deg(theta_z) ;
       
end

    