% SS_PROJ Projects all cell albedo contributions onto the solar cell normal
% in ECEF. This value is equivalent to the total perpendicular irradiance
% reaching the solar cell.
%
% P = ss_proj(re,n,a)
%
% where a is an albedo matrix, n is the cell normal in ECEF frame, and
% re is the nadir in ECEF frame.
%
% $Id: ss_proj.m,v 1.4 2006/05/17 14:39:17 danji Exp $

function P = ss_proj(re,n,a);

CONST.EMR = 6371.01e3;
CONST.AM0 = 1366.9;

P = 0;
grid = zeros(3,1);

[sy sx] = size(a);

for i=1:sy
  for j=1:sx
    if a(i,j) > 0
      % Grid vector in ECEF
      [grid_theta grid_phi] = idx2rad(i,j,sy,sx);
      [grid(1) grid(2) grid(3)] = sph2cart(grid_theta,pi/2-grid_phi,CONST.EMR);
      % Cosine of the angle from solar cell normal to grid LOS vector (re+grid).
      cosphi = n'*(re+grid)/norm(re+grid);
      if cosphi > 0
      angle(i,j) = acos(cosphi);
        P = P + a(i,j)*cosphi;
      end
    end
  end
end

return
