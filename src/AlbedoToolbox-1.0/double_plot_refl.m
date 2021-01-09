% DOUBLE_PLOT_REFL Plot REFL data on subplot with dual view.
%
% double_plot_refl(refl)
%
% refl is either a REFL struct or REFL matrix.
%
% $Id: double_plot_refl.m,v 1.6 2006/05/17 11:13:22 danji Exp $

function double_plot_refl(refl);

subplot(1,2,1);
plot_refl(refl);
subplot(1,2,2);
n = size(refl.data,1);
[sy sx] = size(refl.data);
dy = 180/sy;
lat = [-90+dy/2:dy:90-dy/2]';
plot(lat,mean(refl.data'),'b-',lat,max(refl.data'),'b:',lat,min(refl.data'),'b:');
title('TOMS Reflectivity Data');
xlabel('Latitude [deg]');
ylabel('Reflectivity [E.U.]');
axis([-90+dy 90-dy 0 1]);
legend('Mean','Min/Max');


if nargout > 0
    % Return handle if requested
    h = gcf;
end

return
