function eclipse = calculate_eclipse(xsat_eci,sun_pos_eci)
%% Calculate eclipse
R_EARTH = 6371;
R_SUN=696000;
AU = 149600000;

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
alpha1;
alpha2;
alpha;
eclipse;
end