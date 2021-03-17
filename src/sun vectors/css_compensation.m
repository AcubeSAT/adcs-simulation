%% CSS Compensation
%% Inputs:
% Sun in the ECI frame (sun_eci), the quaternion from ECI to Body
% frame,(q_eci,body), the position of the Satellite in the ECI frame and
% the percentage of the sunlight diffused from earth that interacts with
% the satellite (albedo_perc).
%% Outputs:
% The total sun vector in the Body Frame (total_sun_vector)
%%
% This function is used in the Measurement Function of the Kalman Filter.
% It's purpose is to provide what the measurment of the CSS will be. Hence,
% it follows the same logic with the function CSS Noise, but without adding
% the poisson noise, which of course is random and cannot be modeled.
% Subsequently, the Measurement Function can define the difference between
% the actual measurement (which comes from CSS Noise) and the expected
% measurement (which comes from this function) and adjust accordingly.



function total_sun_vector = css_compensation(sun_eci,q_eci_body,xsat_eci,albedo_perc,~)
    
    sun_eci = sun_eci/norm(sun_eci);
    temp = quatProd(quatconj(q_eci_body'),quatProd([0;sun_eci],q_eci_body));
    sun_body = temp(2:4);
    
    %frame_of_css_1 = roty(0);
    frame_of_css_2 = roty(90);
    frame_of_css_3 = roty(-90);
    frame_of_css_4 = roty(180);
    frame_of_css_5 = rotz(90);
    frame_of_css_6 = rotz(-90);

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
        if current_sun(i) < 0
            current_sun(i) = 0;
        end
        current_albedo(i) = albedo_perc * dot(final_nadir(:,i),[1;0;0]);
        if current_albedo(i) < 0
            current_albedo(i) = 0;
        end
        total_current(i) = current_sun(i) + current_albedo(i);
    end
  
    total_sun_vector = [total_current(1)-total_current(4), total_current(6)-total_current(5),total_current(2)-total_current(3)];
    total_sun_vector = (total_sun_vector/norm(total_sun_vector))'; % Sun vector in body frame
 
end