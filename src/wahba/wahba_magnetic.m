function q_mag_body = wahba_magnetic(mag_ref, mag_body)

    mag_ref = mag_ref / norm(mag_ref);
    mag_body = mag_body / norm(mag_body);

    if dot(mag_ref, mag_body) > 0.999999
        q_mag_body(1) = 1;
        q_mag_body(2:4) = [0; 0; 0];
    elseif dot(mag_ref, mag_body) < -0.999999
        q_mag_body(1) = 0;
        q_mag_body(2:4) = [1; 0; 0];
    else
        a = cross(mag_body, mag_ref);
        q_mag_body(2:4) = a;
        q_mag_body(1) = 1 + dot(mag_ref, mag_body);
        q_mag_body = q_mag_body / norm(q_mag_body);
    end

end
