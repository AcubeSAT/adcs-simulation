% -----------------------------------------------------------------------------
%
%                              procedure calculate_eclipse
%
%  this procedure calculates the extistence of eclipse at any time during in-orbit
%
%   Inputs        :
%    xsat_eci     - Satellite position in ECI frame
%    sun_pos_eci  - Sun position in ECI frame
%
%   Outputs       :
%    eclipse      - Existence (or not) of elcipse in umbral or penumbral conditions
%                   0 -> no eclipse , 1 -> penumbral , 2 -> umbral
%  ----------------------------------------------------------------------------*/

function eclipse = calculate_eclipse(xsat_eci,sun_pos_eci)
%% -------- Constants ------
R_EARTH = 6371;
R_SUN=696000;
AU = 149600000;

%% -------- Algorithm ------
for n=1:size(xsat_eci,2)
    
    % from the penumbral cone geometry
    x1 = R_EARTH*AU/(R_SUN+R_EARTH);
    alpha1(:,n) = pi - acos(R_EARTH/x1)-acos(R_EARTH/norm(xsat_eci(:,n)));

    % from the umbral cone geometry
    x2 = R_EARTH*AU/(R_SUN-R_EARTH);
    alpha2(:,n) = acos(R_EARTH/x2)-acos(R_EARTH/norm(xsat_eci(:,n)));


    alpha(:,n) = pi - acos(sun_pos_eci(:,n)'*xsat_eci(:,n)/(norm(sun_pos_eci(:,n))*norm(xsat_eci(:,n))));


    if alpha2(:,n) < alpha(:,n) & alpha(:,n) < alpha1(:,n)
        % penumbral 
        eclipse(n) = 1;
    elseif alpha(:,n) < alpha2(:,n)
        % umbral
        eclipse(n) = 2;
    else 
        % no eclipse
        eclipse(n) =0;
    end
end
end
