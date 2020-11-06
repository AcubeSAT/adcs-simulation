dt = 1;
npts = 2*3600;
satrec = orbit_init();
[xsat_ecf, vsat_ecf,xsat_eci,vsat_eci, sat_llh,eclipse, mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit,...
    sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,satrec,argpm,nodem,inclm,mm,xnode,xinc] = orbit_sgp4(satrec,dt,npts);
% [xsat_ecf, vsat_ecf,xsat_eci,vsat_eci, sat_llh,eclipse, mag_field_ned,mag_field_eci,mag_field_ecef,mag_field_orbit, sun_pos_ned,sun_pos_eci,sun_pos_ecef,sun_pos_orbit,satrec,argpm,nodem,inclm,mm,xnode,xinc] = orbit_sgp4_offset(satrec,dt,dt,1);

sat_llh(1,:)=sat_llh(1,:)*180/pi;
sat_llh(2,:)=sat_llh(2,:)*180/pi;

t=[0:(npts-1)/dt]*dt;


figure('Position',[0 0 640 1080])
for i=1:3
    subplot(3,1,i);
    hold on;
    plot(t,xsat_eci(i,:), 'LineWidth',2.0, 'Color','blue');
%     legend(['metres'], 'interpreter','latex', 'fontsize',15);
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    ylabel(['$Xeci_' num2str(i) '$'], 'interpreter','latex', 'fontsize',14);
    if (i==1), title('Satellite Position ECI', 'interpreter','latex', 'fontsize',17);end
    hold off;
end
figure('Position',[0 0 640 1080])
subplot(4,1,1);
plot(t,argpm(1,:), 'LineWidth',2.0, 'Color','blue');
ylabel(['$Arg of Perigee$'], 'interpreter','latex', 'fontsize',14);
xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
subplot(4,1,2);
plot(t,nodem(1,:), 'LineWidth',2.0, 'Color','blue');
ylabel(['$RAAN$'], 'interpreter','latex', 'fontsize',14);
xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
subplot(4,1,3);
plot(t,inclm(1,:), 'LineWidth',2.0, 'Color','blue');
ylabel(['$Inclination$'], 'interpreter','latex', 'fontsize',14);
xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
subplot(4,1,4);
plot(t,mm(1,:), 'LineWidth',2.0, 'Color','blue');
ylabel(['$Mean Anomaly$'], 'interpreter','latex', 'fontsize',14);
xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
% figure('Position',[0 0 640 1080])
% for i=1:3
%     subplot(3,1,i);
%     hold on;
%     plot(t,xsat_ecf(i,:), 'LineWidth',2.0, 'Color','blue');
%     legend(['metres'], 'interpreter','latex', 'fontsize',15);
%     xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
%     ylabel(['$Xeci_' num2str(i) '$'], 'interpreter','latex', 'fontsize',14);
%     if (i==1), title('Satellite Position ECEF', 'interpreter','latex', 'fontsize',17);end
%     hold off;
% end

figure('Position',[640 0 640 1080]);
for i=1:3
    subplot(3,1,i);
    hold on;
    plot(t,sat_llh(i,:), 'LineWidth',2.0, 'Color','blue');
    if i~=3
        legend(['degrees'], 'interpreter','latex', 'fontsize',15);
    else
        legend(['metres'], 'interpreter','latex', 'fontsize',15);
    end
    if (i==1)
        ylabel(['latitude'], 'interpreter','latex', 'fontsize',17);
    elseif (i==2)
        ylabel(['longitude'], 'interpreter','latex', 'fontsize',17);
    else
        ylabel(['altitude'], 'interpreter','latex', 'fontsize',17);
    end
    
    if (i==1), title('Latitude, Longitude and Altitude', 'interpreter','latex', 'fontsize',17);end
    hold off;
end


% figure();
% plot(t,sat_llh(1,:)*180/pi);
% 
% figure();
% plot(t,sat_llh(2,:)*180/pi);
% 
% figure();
% plot(t,sat_llh(3,:));
figure('Position',[1280 0 640 1080])
for i=1:3
    subplot(3,1,i);
    hold on;
    plot(t,mag_field_ned(i,:), 'LineWidth',2.0, 'Color','blue');
    legend(['nanoTesla'], 'interpreter','latex', 'fontsize',15);
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    ylabel(['$Bned_' num2str(i) '$'], 'interpreter','latex', 'fontsize',14);
    if (i==1), title('Magnetic field NED', 'interpreter','latex', 'fontsize',17);end
    hold off;
end



figure('Position',[0 0 640 1080])
for i=1:3
    subplot(3,1,i);
    hold on;
    plot(t,sun_pos_eci(i,:), 'LineWidth',2.0, 'Color','blue');
    legend(['au'], 'interpreter','latex', 'fontsize',15);
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    ylabel(['$SUNeci_' num2str(i) '$'], 'interpreter','latex', 'fontsize',14);
    if (i==1), title('Sun Position ECI', 'interpreter','latex', 'fontsize',17);end
end

figure('Position',[0 0 640 360])
for i=1:1
    subplot(1,1,i);
    hold on;
    plot(t,eclipse, 'LineWidth',2.0, 'Color','blue');
%     legend(['au'], 'interpreter','latex', 'fontsize',15);
    xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
    ylabel(['Eclipse'], 'interpreter','latex', 'fontsize',14);
    if (i==1), title('Umbral, Penumbral or no Eclipse', 'interpreter','latex', 'fontsize',17);end
end

% umbral=zeros(length(eclipse));
% penumbral=zeros(length(eclipse));
% normal=zeros(length(eclipse));

% for j=1:length(eclipse)
%     if eclipse(j) == 0
%         normal(j) = 1;
%     elseif eclipse(j) == 1
%         penumbral(j) = 1;
%     else
%         umbral(j) = 1;
%     end
% end
% 
% 
% figure('Position',[0 0 640 1080])
% for i=1:3
%     subplot(3,1,i);
%     hold on;
%     if i==1
%         plot(t,umbral, 'LineWidth',2.0, 'Color','[.61 .51 .74]');
%         ylabel(['Umbral'], 'interpreter','latex', 'fontsize',14);
%     elseif i==2
%         plot(t,penumbral, 'LineWidth',2.0, 'Color','[0, 0.4470, 0.7410]	');
%         ylabel(['Penumbral'], 'interpreter','latex', 'fontsize',14);
%     else
%         plot(t,normal, 'LineWidth',2.0, 'Color','[0.8500, 0.3250, 0.0980]');
%         ylabel(['No Eclipse'], 'interpreter','latex', 'fontsize',14);
%     end
%     %     plot(t,eclipse, 'LineWidth',2.0, 'Color','[0.8500, 0.3250, 0.0980]');
% %     if (i==1), legend({'umbral','penumbral','no eclipse'}, 'interpreter','latex', 'fontsize',15);end
%     xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
%     
% 
%     if (i==1), title('Eclipse', 'interpreter','latex', 'fontsize',17);end
%     hold off;
% end

% figure('Position',[0 0 640 1080])
% for i=1:4
%     subplot(4,1,i);
%     hold on;
%     if i~=4
%         plot(t,sun_pos_eci(i,:), 'LineWidth',2.0, 'Color','blue');
%         legend(['au'], 'interpreter','latex', 'fontsize',15);
%         xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
%         ylabel(['$SUNeci_' num2str(i) '$'], 'interpreter','latex', 'fontsize',14);
%     else
% %         plot(t,umbral, 'LineWidth',2.0, 'Color','[.61 .51 .74]');
% %         plot(t,penumbral, 'LineWidth',2.0, 'Color','[0, 0.4470, 0.7410]	');
% %         plot(t,normal, 'LineWidth',2.0, 'Color','[0.8500, 0.3250, 0.0980]');
%         plot(t,eclipse, 'LineWidth',2.0, 'Color','[0.8500, 0.3250, 0.0980]');
% %         legend(['umbral'],['penumbral'],['no eclipse'], 'interpreter','latex', 'fontsize',15);
%         xlabel('Time [$s$]', 'interpreter','latex', 'fontsize',12);
%         ylabel(['Eclipse (no eclipse, Penumbral, Umbral)'], 'interpreter','latex', 'fontsize',14);
%     end
%     if (i==1), title('Sun Position ECI and Eclipse', 'interpreter','latex', 'fontsize',17);end
%     hold off;
% end

% figure();
% plot(t,mag_field_ned(1,:))
% 
% figure();
% plot(t,mag_field_ned(2,:))
% 
% figure();
% plot(t,mag_field_ned(3,:))

% figure();
% plot(t,mag_field_eci(1,:))
% 
% figure();
% plot(t,mag_field_eci(2,:))
% 
% figure();
% plot(t,mag_field_eci(3,:))

% figure();
% plot(t,mag_field_ecef(1,:))
% 
% figure();
% plot(t,mag_field_ecef(2,:))
% 
% figure();
% plot(t,mag_field_ecef(3,:))

% figure();
% plot(t,mag_field_orbit(1,:))
% 
% figure();
% plot(t,mag_field_orbit(2,:))
% 
% figure();
% plot(t,mag_field_orbit(3,:))

% figure();
% plot(t,sun_pos_ned(1,:))
% 
% figure();
% plot(t,sun_pos_ned(2,:))
% 
% figure();
% plot(t,sun_pos_ned(3,:))

% figure();
% plot(t,sun_pos_orbit(1,:))
% 
% figure();
% plot(t,sun_pos_orbit(2,:))
% 
% figure();
% plot(t,sun_pos_orbit(3,:))

