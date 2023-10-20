%% project: GNSSR_coverage_simulator
% Santiago Ozafrain, Jun 2022. Updated: Oct 2023.
% UIDET-SENyT, Facultad de Ingenieria, UNLP.
%
% Main script - Calculates the position of the specular reflection points
% for a given LEO satellite trajectory. Computes the number of reflections
% captured using an elliptical cone as model for the antenna beam. It also
% calculates relevant parameters of those reflections: incidence angle,
% path loss and Rx antenna gain.

clear;
close all;clc

config

%% Computation of the specular reflection points positions for every GPS
%  satellite in each instant of the time interval.


if load_sp == 0 | ~isfile(sp_filename)
    tic

    wb = waitbar(0, '   Calculating specular reflection points position: 0% completed    ');
    for nn = 1:N
        [satPos, satVel, satID] = gnssconstellation(t(nn), gpsData, GNSSFileType="RINEX");
        gps_pos(nn, :, :) = satPos;
        for SVID = 1:num_gps_sats
            sp_pos_aux(nn,SVID,:) = SpecularReflectionPoint(satPos(SVID,:), leo_pos(nn,:)); % Specular reflection point position
        end

        % progression bar
        time_past = toc;
        time_per_loop = time_past/nn;
        time_left = (N-nn)*time_per_loop;
        time_left_h = floor(time_left/3600);
        time_left_m = floor((time_left - time_left_h*60)/60);
        time_left_s = floor(time_left - time_left_h*3600 - time_left_m*60);
        msg = sprintf('   Calculating specular reflection points position: %i%% completed   \n%i:%i:%i remaining', floor(nn/N*100), time_left_h, time_left_m, time_left_s);
        waitbar(nn/N, wb, msg)
    end
    close(wb)
    save(sp_filename)
else
    load(sp_filename)
end

%% Computation of reflections within the antenna footprint
% It checks which SP are located inside the antenna footprint. That
% information is stored as a flag in sp_stat for each instant and SVID,
% with 1 indicating that that reflection is being captured by the LEO
% satellite.

[X, Y, Z] = gen_beam(theta_x, theta_y, 1e5, 8e5); % Antenna beam generation

tic
wb = waitbar(0, '   Counting specular reflection points within antenna footprint: 0% completed   ');

for nn = 10 : N % loops for the total time interval

    leo_vel = [mean(diff(leo_pos(nn-9:nn,1))), mean(diff(leo_pos(nn-9:nn,2))), mean(diff(leo_pos(nn-9:nn,3)))]'; % LEO satellite velocity vector

    % Unit vectors to rotate and translate ECEF to RTS
    wp = leo_pos(nn,:)/norm(leo_pos(nn,:)); vp = cross(wp,leo_vel); vp = vp/norm(vp); up = cross(vp,wp); up = up/norm(up);
    upp = up*cos(theta_x0*pi/180) - wp*sin(theta_x0*pi/180);
    vpp = -up*sin(theta_x0*pi/180)*sin(theta_y0*pi/180) + vp*cos(theta_y0*pi/180) - wp*cos(theta_x0*pi/180)*sin(theta_y0*pi/180);
    wpp = up*sin(theta_x0*pi/180)*cos(theta_y0*pi/180) + vp*sin(theta_y0*pi/180) + wp*cos(theta_x0*pi/180)*cos(theta_y0*pi/180);


    for SVID = 1:num_gps_sats % loops for every GPS satellite

        % convert SP location to the coordinate system fixed in the LEO
        % satellite (RTS)
        sp_pos_aux_rts(nn-9,SVID,:) = [(squeeze(sp_pos_aux(nn,SVID,:))'-leo_pos(nn,:))*upp' (squeeze(sp_pos_aux(nn,SVID,:))'-leo_pos(nn,:))*vpp' (squeeze(sp_pos_aux(nn,SVID,:))'-leo_pos(nn,:))*wpp']';

        if isinfp(sp_pos_aux_rts(nn-9,SVID,:), X, Y, Z) % checks if SP is within antenna footprint

            % reflection captured flag
            sp_stat(nn-9,SVID) = 1;

            % incidence angle of the captured reflection
            sp_inc_ang_aux(nn-9, SVID) = acos(-wp*(squeeze(sp_pos_aux(nn,SVID,:))-leo_pos(nn,:)')/norm(squeeze(sp_pos_aux(nn,SVID,:))'-leo_pos(nn,:)));


            % Rx antena gain for the captured SP location
            el_sp(nn-9, SVID) = acos(-[0 0 1]*squeeze(sp_pos_aux_rts(nn-9,SVID,:))/norm(squeeze(sp_pos_aux_rts(nn-9,SVID,:))));
            az_sp_aux = atan2(sp_pos_aux_rts(nn-9,SVID,2), sp_pos_aux_rts(nn-9,SVID,1));
            az_sp_aux(az_sp_aux<0) = az_sp_aux + 2*pi;
            az_sp(nn-9, SVID) = az_sp_aux;

            [~, el_idx] = min(abs(el-el_sp(nn-9, SVID)*180/pi)); % Find index to gain value closer to SP el angle
            [~, az_idx] = min(abs(az-az_sp(nn-9, SVID)*180/pi)); % Find index to gain value closer to SP az angle
            G_ant_aux(nn-9, SVID) = G0(el_idx,az_idx);

            % path loss of the captured reflection (non-coherent/diffuse reflection)
            sp_path_loss_aux(nn-9, SVID) = 10*log10(norm(squeeze(sp_pos_aux(nn,SVID,:))'-leo_pos(nn,:))^2*norm(squeeze(sp_pos_aux(nn,SVID,:)-gps_pos(nn,SVID,:)))^2);
        end
    end

    % progression bar
    time_past = toc;
    time_per_loop = time_past/nn;
    time_left = (N-nn)*time_per_loop;
    time_left_h = floor(time_left/3600);
    time_left_m = floor((time_left - time_left_h*60)/60);
    time_left_s = floor(time_left - time_left_h*3600 - time_left_m*60);
    msg = sprintf('   Counting specular reflection points within antenna footprint: %i%% completed\n   %i:%i:%i remaining', floor(nn/N*100), time_left_h, time_left_m, time_left_s);
    waitbar(nn/N, wb, msg)
end
close(wb)

% Cell objects store values only for captured reflections
for SVID = 1:num_gps_sats
    sp_pos{SVID} = sp_pos_aux(sp_stat(:,SVID)==1,SVID,:);
    sp_pos_rts{SVID} = sp_pos_aux_rts(sp_stat(:,SVID)==1,SVID,:);
    sp_inc_ang{SVID} = sp_inc_ang_aux(sp_stat(:,SVID)==1,SVID);
    sp_path_loss{SVID} = sp_path_loss_aux(sp_stat(:,SVID)==1,SVID);
    G_ant{SVID} = G_ant_aux(sp_stat(:,SVID)==1,SVID);
end

sp_num = sum(sp_stat,2); % Number of reflections captured in each instant
clear sp_inc_ang_aux sp_path_loss_aux G_ant_aux % Clear auxiliar variables

%% Analysis with limited time registers

sp_pres = max(sp_stat,[],2); % flag: 1 if there is at least 1 reflection captured

for nn = 1:length(1:length(sp_pres)-reg_len_samples)
    % sums the time that there is at least 1 reflection captured in the
    % register of reg_len_samples duration (reg_len in seconds)
    sp_time(nn) = sum(sp_pres(nn:nn+reg_len_samples))*dt;
end

%% Reflections over Argentina

for SVID = 1:num_gps_sats
    inArg_cell = zeros(size(squeeze(sp_pos{SVID})));
    inArg_array = zeros(size(sp_stat(sp_stat(:,SVID)==1,SVID)));

    if any(squeeze(sp_pos{SVID}))
        inArg_cell = isinArg(squeeze(sp_pos{SVID})); % Checks if the reflection is over Argentina
        inArg_array = isinArg(squeeze(sp_pos_aux(sp_stat(:,SVID)==1,SVID,:)));
    end
    sp_pos_arg{SVID} = squeeze(sp_pos{SVID}).*inArg_cell;
    sp_stat_arg(sp_stat(:,SVID)==1,SVID) = sp_stat(sp_stat(:,SVID)==1,SVID).*inArg_array;
end

sp_num_arg = sum(sp_stat_arg,2); % Number of reflections captured over Argentina in each instant


%% Plot results
plot_results

%%
% save('files/res00.mat')