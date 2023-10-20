%% project: GNSSR_coverage_simulator
%
%
% Script to plot simulation results.


close all


% Ellipsoidal Earth model WGS48
E=referenceEllipsoid(7030);
xr=E.SemimajorAxis;
yr=E.SemimajorAxis;
zr=E.SemiminorAxis;
[xe,ye,ze] = ellipsoid(0,0,0,xr,yr,zr,750);

% GPS sat trajectory
SVID = 1; % GPS sat ID number to plot

% LEO sat position and antenna footprint
t_snapshot = 100; % instant selected to plot
leo_vel = [mean(diff(leo_pos(t_snapshot-9:t_snapshot, 1))), mean(diff(leo_pos(t_snapshot-9:t_snapshot, 2))), mean(diff(leo_pos(t_snapshot-9:t_snapshot, 3)))]';

% Antenna beam generation
[X, Y, Z] = gen_beam(theta_x, theta_y, 0, leo_alt*1.1);

% Antenna beam rotation to plot in ECEF
u = [1 0 0]'; v = [0 1 0]'; w = [0 0 1]';
wp = -leo_pos(t_snapshot, :)'/norm(leo_pos(t_snapshot, :)); vp = cross(leo_vel,wp); vp = vp/norm(vp); up = cross(vp,wp); up = up/norm(up);
R = [u'*up v'*up w'*up;
    u'*vp v'*vp w'*vp;
    u'*wp v'*wp w'*wp];
R = Rt * R;
Xr = R(1,1)*X + R(2,1)*Y - R(3,1)*Z + leo_pos(t_snapshot, 1);
Yr = R(1,2)*X + R(2,2)*Y - R(3,2)*Z + leo_pos(t_snapshot, 2);
Zr = R(1,3)*X + R(2,3)*Y - R(3,3)*Z + leo_pos(t_snapshot, 3);

%% plot: Bistatic radar geometry: GPS satellite orbit, LEO position and
% antenna beam, specular reflection point position and vector pointing in 
% the moving direction.

figure(1);
plot3(gps_pos(:,SVID,1)*1e-3,gps_pos(:,SVID,2)*1e-3,gps_pos(:,SVID,3)*1e-3);hold on;
plot3(gps_pos(t_snapshot,SVID,1)*1e-3,gps_pos(t_snapshot,SVID,2)*1e-3,gps_pos(t_snapshot,SVID,3)*1e-3,'*r');
text(gps_pos(t_snapshot,SVID,1)*1e-3+1e3,gps_pos(t_snapshot,SVID,2)*1e-3+1e3,gps_pos(t_snapshot,SVID,3)*1e-3+1e3,sprintf('PRN: %i', SVID));
plot3(sp_pos_aux(t_snapshot,SVID,1)*1e-3,sp_pos_aux(t_snapshot,SVID,2)*1e-3,sp_pos_aux(t_snapshot,SVID,3)*1e-3,'*r');
plot3(sp_pos_aux(t_snapshot,SVID,1)*1e-3*sp_stat(t_snapshot,SVID),sp_pos_aux(t_snapshot,SVID,2)*sp_stat(t_snapshot,SVID)*1e-3,sp_pos_aux(t_snapshot,SVID,3)*sp_stat(t_snapshot,SVID)*1e-3,'*g');
text(sp_pos_aux(t_snapshot,SVID,1)*1e-3+1e2,sp_pos_aux(t_snapshot,SVID,2)*1e-3+1e2,sp_pos_aux(t_snapshot,SVID,3)*1e-3+1e2,'SP','FontSize',10);
surf(xe*1e-3,ye*1e-3,ze*1e-3,'FaceAlpha',0.75,'EdgeColor','none');axis equal;
surf(Xr*1e-3,Yr*1e-3,Zr*1e-3,'EdgeColor','none'); %
line([0 leo_pos(t_snapshot, 1)*1e-3], [0 leo_pos(t_snapshot, 2)*1e-3], [0 leo_pos(t_snapshot, 3)*1e-3])
text(leo_pos(t_snapshot,1)*1e-3+1e2,leo_pos(t_snapshot,2)*1e-3+1e2,leo_pos(t_snapshot,3)*1e-3+1e2,'LEO','FontSize',10);

leo_vel=leo_vel/norm(leo_vel)*1e6;
quiver3(leo_pos(t_snapshot, 1)*1e-3,leo_pos(t_snapshot, 2)*1e-3,leo_pos(t_snapshot, 3)*1e-3,leo_vel(1)*1e-3,leo_vel(2)*1e-3,leo_vel(3)*1e-3,'color','r','LineWidth',1.2);
hold off
xlabel('x [km]');ylabel('y [km]');zlabel('z [km]');grid on
title('Bistatic radar scenario')

%% plot: Reflections captured in the time interval
figure(2);
surf(xe*1e-3,ye*1e-3,ze*1e-3,'EdgeColor','none');axis equal;grid on
hold on
for SVID = 1:num_gps_sats
    plot3(sp_pos{SVID}(:, 1)*1e-3,sp_pos{SVID}(:, 2)*1e-3,sp_pos{SVID}(:, 3)*1e-3,'.g')
end
xlabel('x [km]');ylabel('y [km]');zlabel('z [km]');grid on
title('Specular reflection points of captured signals - ECEF')
hold off

%% plot: Specular reflection points position RTS

figure(3);
for SVID = 1:num_gps_sats
    plot3(sp_pos_aux_rts(:,SVID,1)*1e-3, sp_pos_aux_rts(:,SVID,2)*1e-3, sp_pos_aux_rts(:,SVID,3)*1e-3,'.r');hold on
    plot3(sp_pos_rts{SVID}(:, 1)*1e-3, sp_pos_rts{SVID}(:, 2)*1e-3, sp_pos_rts{SVID}(:, 3)*1e-3,'.g');
end
surf(X*1e-3,Y*1e-3,Z*1e-3,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('x [km]');ylabel('y [km]');zlabel('z [km]');grid on
title('Specular reflection points position relative to LEO satellite')
hold off
%% plot: Number of reflections captured simultaneously
figure(4);

plot(t0/3600, sp_num);grid on
title('Number of reflected signals captured')
axis([t0(1)/3600 t0(end)/3600 -0.2 max(sp_num)+0.2])
xlabel('time [hs]')
%% plot: Histogram of incindence angle of captured reflections
figure(5);
angs = [];
for SVID = 1:num_gps_sats
    angs = [angs; sp_inc_ang{SVID}];
end
angs=angs*180/pi;
histogram(angs); xlabel('\theta_{inc} [degs]');grid on
title('Histogram - Incidence angle of captured reflected signals')

%% plot: Histogram of path loss of captured reflections
figure(6);
path_loss = [];
for SVID = 1:num_gps_sats
    path_loss = [path_loss; sp_path_loss{SVID}];
end
histogram(path_loss); xlabel('L_r [dB]');grid on
title("Histogram - Path loss of captured reflected signals - mean value = "+num2str(mean(path_loss)) + " dB")

%% plot: Histogram of receiver antenna gain of captured reflections
figure(7);
G = [];
for SVID = 1:num_gps_sats
    G = [G; G_ant{SVID}];
end
histogram(G); xlabel('G_{Rx} [dB]');grid on
title("Histogram - Receiver antenna gain for each SP - mean value = "+num2str(mean(G)) + " dB")

%% plot: Histogram of normalized power level of captured reflections
figure(8);
histogram(G-path_loss); xlabel('G_{Rx} - L_r [dB]');grid on
title("Histogram - Normalized power level - mean value = "+num2str(mean(G-path_loss)) + " dB")

%% plot: Histogram of time with reflections present in the register
figure(9);
histogram(sp_time);xlabel('time with reflections present [secs]');grid on
title(sprintf('Histogram - Time with at least 1 reflection in a %i seconds register', reg_len))

%% plot: Number of reflections captured simultaneously over Argentina
figure(10);
plot(t0/3600, sp_num_arg);grid on
title('Number of reflected signals captured over Argentina')
axis([t0(1)/3600 t0(end)/3600 -0.2 max(sp_num)+0.2])
xlabel('time [hs]')
%% plot: Reflections captured in the time interval over Argentina
figure(11);
surf(xe,ye,ze,'EdgeColor','none');axis equal;grid on
hold on
for SVID = 1:num_gps_sats
    plot3(sp_pos_arg{SVID}(:, 1),sp_pos_arg{SVID}(:, 2),sp_pos_arg{SVID}(:, 3),'.g')
end
xlabel('x');ylabel('y');zlabel('z');
hold off
title('Specular reflection points of captured signals over Argentina - ECEF')

%% plot: Earth globe showing reflections over Argentina
uif = uifigure;
g = geoglobe(uif,"Basemap","landcover");
lla = [];
for SVID = 1:num_gps_sats
    lla = [lla; ecef2lla(sp_pos_arg{SVID})];
end
p=geoplot3(g,lla(:,1),lla(:,2),lla(:,3),'c','LineWidth',2);

g.Terrain = 'none';
p.HeightReference = 'ellipsoid';
