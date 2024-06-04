%% Parameter initialization
% -----------------------------------------------------------------------------
%     dt      - timestep for the orbit propagator
%     npts    - total simulation seconds
%  ----------------------------------------------------------------------------*/

dt = 1;
npts = 3 * 3600;

%% Orbit propagator
satrec = orbit_init();
[xsat_ecf, vsat_ecf, xsat_eci, vsat_eci, sat_llh, eclipse, mag_field_ned, mag_field_eci, mag_field_ecef, mag_field_orbit, ...
    sun_pos_ned, sun_pos_eci, sun_pos_ecef, sun_pos_orbit, satrec, argpm, nodem, inclm, mm, xnode, xinc] = orbit_sgp4(satrec, dt, npts);
%[xsat_ecf, vsat_ecf,xsat_eci,vsat_eci, sat_llh,eclipse, mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit, ...
%sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,satrec,argpm,nodem,inclm,mm,xnode,xinc] = orbit_sgp4_offset(satrec,dt,npts,1000);

% Deg2Rad
sat_llh(1, :) = sat_llh(1, :) * 180 / pi;
sat_llh(2, :) = sat_llh(2, :) * 180 / pi;

t = [0:(npts - 1) / dt] * dt;

%% Plots

%% Satellite Position in ECI frame

% figure('Position',[0 0 640 1080])
figure()
for i = 1:3
    subplot(3, 1, i);
    hold on;
    plot(t, xsat_eci(i, :), 'LineWidth', 2.0, 'Color', 'blue');
    xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
    ylabel(['$Xeci_', num2str(i), '$'], 'interpreter', 'latex', 'fontsize', 14);
    if (i == 1), title('Satellite Position ECI', 'interpreter', 'latex', 'fontsize', 17); end
    hold off;
end

%% Orbital Parameters

% figure('Position',[0 0 640 1080])
figure()
subplot(4, 1, 1);
plot(t, argpm(1, :), 'LineWidth', 2.0, 'Color', 'blue');
ylabel(['$Arg of Perigee$'], 'interpreter', 'latex', 'fontsize', 14);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
subplot(4, 1, 2);
plot(t, nodem(1, :), 'LineWidth', 2.0, 'Color', 'blue');
ylabel(['$RAAN$'], 'interpreter', 'latex', 'fontsize', 14);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
subplot(4, 1, 3);
plot(t, inclm(1, :), 'LineWidth', 2.0, 'Color', 'blue');
ylabel(['$Inclination$'], 'interpreter', 'latex', 'fontsize', 14);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
subplot(4, 1, 4);
plot(t, mm(1, :), 'LineWidth', 2.0, 'Color', 'blue');
ylabel(['$Mean Anomaly$'], 'interpreter', 'latex', 'fontsize', 14);
xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);

%% Satellite Position in Latitude, Longitude, Altitude

% figure('Position',[640 0 640 1080]);
figure()
for i = 1:3
    subplot(3, 1, i);
    hold on;
    plot(t, sat_llh(i, :), 'LineWidth', 2.0, 'Color', 'blue');
    if i ~= 3
        legend(['degrees'], 'interpreter', 'latex', 'fontsize', 15);
    else
        legend(['metres'], 'interpreter', 'latex', 'fontsize', 15);
    end
    if (i == 1)
        ylabel(['latitude'], 'interpreter', 'latex', 'fontsize', 17);
    elseif (i == 2)
        ylabel(['longitude'], 'interpreter', 'latex', 'fontsize', 17);
    else
        ylabel(['altitude'], 'interpreter', 'latex', 'fontsize', 17);
    end

    if (i == 1), title('Latitude, Longitude and Altitude', 'interpreter', 'latex', 'fontsize', 17); end
    hold off;
end

%% Magnetic Field in NED frame

% figure('Position',[1280 0 640 1080])
figure()
for i = 1:3
    subplot(3, 1, i);
    hold on;
    plot(t, mag_field_ned(i, :), 'LineWidth', 2.0, 'Color', 'blue');
    legend(['nanoTesla'], 'interpreter', 'latex', 'fontsize', 15);
    xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
    ylabel(['$Bned_', num2str(i), '$'], 'interpreter', 'latex', 'fontsize', 14);
    if (i == 1), title('Magnetic field NED', 'interpreter', 'latex', 'fontsize', 17); end
    hold off;
end

%% Sun Position in ECI frame

% figure('Position',[0 0 640 1080])
figure()
for i = 1:3
    subplot(3, 1, i);
    hold on;
    plot(t, sun_pos_eci(i, :), 'LineWidth', 2.0, 'Color', 'blue');
    legend(['au'], 'interpreter', 'latex', 'fontsize', 15);
    xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
    ylabel(['$SUNeci_', num2str(i), '$'], 'interpreter', 'latex', 'fontsize', 14);
    if (i == 1), title('Sun Position ECI', 'interpreter', 'latex', 'fontsize', 17); end
end

%% Eclipse calculation

% figure('Position',[0 0 640 360])
figure()
for i = 1:1
    subplot(1, 1, i);
    hold on;
    plot(t, eclipse, 'LineWidth', 2.0, 'Color', 'blue');
    %     legend(['au'], 'interpreter','latex', 'fontsize',15);
    xlabel('Time [$s$]', 'interpreter', 'latex', 'fontsize', 12);
    ylabel(['Eclipse'], 'interpreter', 'latex', 'fontsize', 14);
    if (i == 1), title('Umbral [2], Penumbral [1] or no Eclipse [0]', 'interpreter', 'latex', 'fontsize', 17); end
end
