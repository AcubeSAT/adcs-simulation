% GRIDANGLE Calculate angle between two grid index pairs.
%
% rho = gridangle(i1,j1,i2,j2,sy,sx)
%
% $Id: gridangle.m,v 1.3 2006/05/17 14:39:17 danji Exp $

function rho = gridangle(i1,j1,i2,j2,sy,sx);

[theta1,phi1] = idx2rad(i1,j1,sy,sx);
[theta2,phi2] = idx2rad(i2,j2,sy,sx);

rho = acos(sin(phi1)*sin(phi2)*cos(theta1-theta2)+cos(phi1)*cos(phi2));

return
