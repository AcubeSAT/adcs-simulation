function total_sun_vector = css_compensation(sun_eci,q_eci_body,xsat_eci,albedo_perc,lambda)
    
    sun_eci = sun_eci/norm(sun_eci);
    temp = quatProd(quatconj(q_eci_body'),quatProd([0;sun_eci],q_eci_body));
    sun_body = temp(2:4);
    
    %f1 = roty(0);
    f2 = roty(90);
    f3 = roty(-90);
    f4 = roty(180);
    f5 = rotz(90);
    f6 = rotz(-90);

    xsat_eci = xsat_eci/norm(xsat_eci);
    xsat_body = quatProd(quatconj(q_eci_body'),quatProd([0;xsat_eci],q_eci_body));
    nadir = -xsat_body(2:4);
    
    %q1 = dcm2quat(f1);
    q2 = dcm2quat(f2);
    q3 = dcm2quat(f3);
    q4 = dcm2quat(f4);
    q5 = dcm2quat(f5);
    q6 = dcm2quat(f6);

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
    end
  
    total_sun_vector = [total_current(1)-total_current(4), total_current(6)-total_current(5),total_current(2)-total_current(3)];
    total_sun_vector = (total_sun_vector/norm(total_sun_vector))'; % Sun vector in body frame
 
end