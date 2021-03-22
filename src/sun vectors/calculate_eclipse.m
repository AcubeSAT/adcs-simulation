% -----------------------------------------------------------------------------
%
%                              procedure calculate_eclipse
%
%  this procedure extracts all required data from the TLE file
%  and initializes the SGP4 orbit propagator
%
%   inputs        :
%   xsat_eci      - Satellite position in ECI frame
%   sun_pos_eci   - Sun position in ECI frame
%
%   outputs       :
%     eclipse     - Ecistence (or not) of elcipse in umbral or penumbral conditions
%  ----------------------------------------------------------------------------*/

function eclipse = calculate_eclipse(xsat_eci,sun_pos_eci)
%% -------- Constants ------
R_EARTH = 6371;
R_SUN=696000;
AU = 149600000;

%% -------- Algorithm ------
for n=1:size(xsat_eci,2)
    
    x1 = R_EARTH*AU/(R_SUN+R_EARTH);
    alpha1(:,n) = pi - acos(R_EARTH/x1)-acos(R_EARTH/norm(xsat_eci(:,n)));

    x2 = R_EARTH*AU/(R_SUN-R_EARTH);
    alpha2(:,n) = acos(R_EARTH/x2)-acos(R_EARTH/norm(xsat_eci(:,n)));

    alpha(:,n) = pi - acos(sun_pos_eci(:,n)'*xsat_eci(:,n)/(norm(sun_pos_eci(:,n))*norm(xsat_eci(:,n))));

    if alpha2(:,n) < alpha(:,n) & alpha(:,n) < alpha1(:,n)
        eclipse(n) = 1;
    elseif alpha(:,n) < alpha2(:,n)
        eclipse(n) = 2;
    else 
        eclipse(n) =0;
    end
end
end
