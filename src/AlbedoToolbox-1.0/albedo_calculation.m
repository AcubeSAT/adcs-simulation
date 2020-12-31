function [total_albedo] = albedo_calculation(i,xsat_ecf,sun_pos_ecef,dt,tf,refl)
    
   % new_refl = replace_nan(refl,'ga050101-051231.mat');

    sun_pos_ecef_m = 1.495978707e+11 * sun_pos_ecef;

    a = albedo(xsat_ecf(:,i),sun_pos_ecef_m(:,i),refl);

    total_albedo = sum(sum(a)); % summing all elements

end

