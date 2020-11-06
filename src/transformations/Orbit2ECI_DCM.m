function R_star = Orbit2ECI_DCM(nodeo,incl,argp_mo)
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

R_star = [-R(1,:); R(3,:); R(2,:)]';

end