function q = q_sun_body(sun_eci, q_eci_body,sun_desired)

    sun_eci = sun_eci / norm(sun_eci);
    q_sun_body = quatProd(quatconj(q_eci_body'), quatProd([0; sun_eci]', q_eci_body'));
    sun_body = q_sun_body(2:4);
    sun_desired = sun_desired / norm(sun_desired);
    

    if dot(sun_body, sun_desired) > 0.999999
        q_sun_body(1) = 1;
        q_sun_body(2:4) = [0; 0; 0];
    elseif dot(sun_body, sun_desired) < -0.999999
        q_sun_body(1) = 0.2;
        q_sun_body(2:4) = [-0.4; -0.4; -0.8];
    else
        a = cross(sun_body, sun_desired);
        q_sun_body(2:4) = a;
        q_sun_body(1) = 1 + dot(sun_body, sun_desired);
        q_sun_body = q_sun_body / norm(q_sun_body);
    end
    q = q_sun_body';

end
