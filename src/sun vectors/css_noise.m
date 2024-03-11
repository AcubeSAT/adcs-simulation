%% CSS Noise Function
%% Inputs: 
% sun_eci = Sun in the ECI frame 
% q_eci_body = the quaternion from ECI to Body frame
% xsat_eci = the position of the Satellite in the ECI frame
% albedo_perc = the percentage of the sunlight diffused from earth that
% interacts with the satellite
% labda = the poisson parameter used to add noise
%% Outputs: 
% total_sun_vector = the total sun vector, i.e. including both sun, albedo, 
% and error measurement in the body frame
%% 
% The function finds where each of the 6 Coarse sun sensors look at the
% giving time by turning the Body frame to an angle according to the
% position of the sensor, and then calculates how much albedo each sensor 
% receives and adds this as a noise in the measurement. All albedo radiation 
% is asumed to come from Nadir. Furthermore, in all 6 measurements, 
% a Poisson noise is added, since CSS are not 100% accurate. 
% The measurements of the CSS are 6 currents, from which we can calculate 
% the total sun vector in the Body Frame. If one of the currents is negative, 
% it means that this sensor does not see the sun.
% For more details, you can refer to DDJF_AOCS file.

function total_sun_vector = css_noise(sun_eci,q_eci_body,xsat_eci,albedo_perc,lambda)
    
    sun_eci = sun_eci/norm(sun_eci);
    temp = quatProd(quatconj(q_eci_body'),quatProd([0;sun_eci],q_eci_body));
    sun_body = temp(2:4);
    
    %frame_of_css_1 = roty(0); the first frame is the same with the Body!
    % frame_of_css_2 = roty(90);
    % frame_of_css_3 = roty(-90);
    % frame_of_css_4 = roty(180);
    % frame_of_css_5 = rotz(90);
    % frame_of_css_6 = rotz(-90);

    frame_of_css_2 = rot_y(90);
    frame_of_css_3 = rot_y(-90);
    frame_of_css_4 = rot_y(180);
    frame_of_css_5 = rot_z(90);
    frame_of_css_6 = rot_z(-90);

    xsat_eci = xsat_eci/norm(xsat_eci);
    xsat_body = quatProd(quatconj(q_eci_body'),quatProd([0;xsat_eci],q_eci_body));
    nadir = -xsat_body(2:4);
    
    %q1 = dcm2quat(f1);
    q2 = dcm2quat(frame_of_css_2);
    q3 = dcm2quat(frame_of_css_3);
    q4 = dcm2quat(frame_of_css_4);
    q5 = dcm2quat(frame_of_css_5);
    q6 = dcm2quat(frame_of_css_6);

    temp_sun(:,1) = [0;sun_body];
    temp_sun(:,2) = quatProd(quatconj(q2),quatProd([0;sun_body],q2));
    temp_sun(:,3) = quatProd(quatconj(q3),quatProd([0;sun_body],q3));
    temp_sun(:,4) = quatProd(quatconj(q4),quatProd([0;sun_body],q4));
    temp_sun(:,5) = quatProd(quatconj(q5),quatProd([0;sun_body],q5));
    temp_sun(:,6) = quatProd(quatconj(q6),quatProd([0;sun_body],q6));

    temp_nadir(:,1) = [0;nadir];
    temp_nadir(:,2) = quatProd(quatconj(q2),quatProd([0;nadir],q2));
    temp_nadir(:,3) = quatProd(quatconj(q3),quatProd([0;nadir],q3));
    temp_nadir(:,4) = quatProd(quatconj(q4),quatProd([0;nadir],q4));
    temp_nadir(:,5) = quatProd(quatconj(q5),quatProd([0;nadir],q5));
    temp_nadir(:,6) = quatProd(quatconj(q6),quatProd([0;nadir],q6));

    final_sun = temp_sun(2:4,:);
    final_nadir = temp_nadir(2:4,:); 
    
    current_sun = zeros(6,1);
    current_albedo = zeros(6,1);
    total_current = zeros(6,1);
    
    for i = 1:6
        current_sun(i) = dot(final_sun(:,i),[1;0;0]);
        sign=randi([0 1]); 
        if sign==0 
           sign=-1; 
        end
        if current_sun(i) < 0
            current_sun(i) = 0;
        end
        current_albedo(i) = albedo_perc * dot(final_nadir(:,i),[1;0;0]);
        if current_albedo(i) < 0
            current_albedo(i) = 0;
        end
        total_current(i) = current_sun(i) + current_albedo(i);
         total_current(i) = total_current(i) + total_current(i)*sign*0.01*poissrnd(lambda);
    end
  
    total_sun_vector = [total_current(1)-total_current(4), total_current(6)-total_current(5),total_current(2)-total_current(3)];
    total_sun_vector = (total_sun_vector/norm(total_sun_vector))'; % Sun vector in body frame
end
