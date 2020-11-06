function ref_vector_orbit = ECI2Orbit(ref_vector_eci,nodeo,incl,argp_mo)
R=zeros(3,3);
R(1,1) = -sin(argp_mo) * sin(nodeo) * cos(incl) + cos(argp_mo) * cos(nodeo);
R(1,2) = sin(argp_mo) * cos(incl) * cos(nodeo) + sin(nodeo) * cos(argp_mo);
R(1,3) = sin(incl) * sin(argp_mo);

R(2,1) = -sin(argp_mo) * cos(nodeo) - sin(nodeo) * cos(incl) * cos(argp_mo);
R(2,2) = -sin(argp_mo) * sin(nodeo) + cos(incl) * cos(argp_mo) * cos(nodeo);
R(2,3) = sin(incl) * cos(argp_mo);

R(3,1) = sin(incl) * sin(nodeo);
R(3,2) = -sin(incl) * cos(nodeo);
R(3,3) = cos(incl);

%X Axis with a minus to point at nadir
for n=1:size(ref_vector_eci,2)
    ref_vector_orbit(1,n)=-(R(1,1)*ref_vector_eci(1,n)+R(1,2)*ref_vector_eci(2,n)+R(1,3)*ref_vector_eci(3,n));
    ref_vector_orbit(2,n)=R(3,1)*ref_vector_eci(1,n)+R(3,2)*ref_vector_eci(2,n)+R(3,3)*ref_vector_eci(3,n);
    ref_vector_orbit(3,n)=R(2,1)*ref_vector_eci(1,n)+R(2,2)*ref_vector_eci(2,n)+R(2,3)*ref_vector_eci(3,n);
end
R;
end