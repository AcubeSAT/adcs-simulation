% EARTH_MEAN_REFL Calculate mean Earth reflectivity of REFL data.
% Each data point is weighed with respect to the cell area of the
% associated cell.
%
% total_refl = earth_mean_refl(refl)
%
% refl is a REFL struct.
%
% $Id: earth_mean_refl.m,v 1.7 2006/05/17 14:39:17 danji Exp $

function total_refl = earth_mean_refl(refl);

CONST.EMR = 6371.01e3;

total_area = 4*pi*CONST.EMR^2;
total_refl = 0;
for i=1:180
  for j=1:288
    weight = cellarea(i,j)/total_area;
    total_refl = total_refl + weight*refl.data(i,j);
  end
end

return
