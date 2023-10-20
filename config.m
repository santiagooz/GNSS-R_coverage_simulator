%% project: GNSSR_coverage_simulator
%
% Configuration script.

%% Input parameteres

time_interval = 1 * 24 * 3600;   % total time interval for the simulation [secs]
dt = 1;                          % time resolution [integer secs]


% ANTENNA BEAM - aperture angles can be extracted from an antenna pattern
% file or set as input parameters
ant_filename = 'files/Realized_Gain_Table_L1.xlsx'; % filename: antenna gain pattern
theta_x0 = 0; % Maximum gain direction along track [degrees]
theta_y0 = 0; % Maximum gain direction across track [degrees]
theta_x = 15; theta_y = 45; % aperture angles in case there is no antenna pattern file present [degrees]


% GPS SATELLITES
rinex_filename = 'files/ADIS00ETH_R_20232470000_01D_GN.rnx'; % filename: GPS rinex file


% LEO SATELLITE
leo_filename = 'files/PosW.42'; % filename: LEO orbit ECEF
leo_alt = 600e3; % nominal LEO sat altitude


% REGISTER LENGTH
reg_len = 60; % register length in [secs]


load_sp = 1; % load SP locations file flag: 0 calculates SP locations - 1 loads file (saves time)
sp_filename = strcat("files/sp_locations_"+"dt_"+num2str(dt)+"_dur_"+num2str(time_interval)+".mat"); % sp locations filename to save/load

%%
% TIME INTERVAL
t_init =  datetime(2023,9,1,14,52,58); % initial time and date
t0 = 1 : dt : time_interval;
t = t_init + t0/24/3600;
N = length(t); % number of instants in the time interval
reg_len_samples = floor(reg_len/dt); % register length in samples

% Rotation matrix for antenna beam direction
Rt = [1 0 0; 0 cos(theta_y0*pi/180) sin(theta_y0*pi/180) ; 0 -sin(theta_y0*pi/180) cos(theta_y0*pi/180)]*[cos(theta_x0*pi/180) 0 sin(theta_x0*pi/180);0 1 0 ; -sin(theta_x0*pi/180) 0 cos(theta_x0*pi/180)];


% GPS satellites orbital parameters
data = rinexread(rinex_filename);
gpsData = data.GPS;
[~,satIdx] = unique(gpsData.SatelliteID);
gpsData = gpsData(satIdx,:);
num_gps_sats = length(satIdx); % number of GPS satellites in the simulation

% LEO satellite trajectory
LEO = load(leo_filename);
leo_pos = LEO(t0,:);

% Antenna pattern
if isfile(ant_filename)
    G0 = xlsread(ant_filename); G0 = 10*log10(G0);
    [el, az, theta_x, theta_y] = antenna_parameters(G0, 2);
else
    G0 = 0; el = 0; az = 0;
end

% output and auxiliar arrays initialization - RTS stands for "relative to
% satellite": reference system fixed to the LEO satellite position
gps_pos = zeros(N, num_gps_sats, 3);          % GPS satellites' position - ECEF
sp_pos_aux = zeros(N, num_gps_sats, 3);       % Specular reflection points' position - ECEF
sp_pos = cell(num_gps_sats,1);                % Specular reflection points' position - ECEF (cell)
sp_pos_arg = cell(num_gps_sats,1);            % Specular reflection points' over Argentina position - ECEF (cell)
sp_pos_aux_rts = zeros(N, num_gps_sats, 3);   % Specular reflection points' position - RTS
sp_pos_rts = cell(num_gps_sats,1);            % Specular reflection points' position - RTS (cell)
sp_stat = zeros(N,num_gps_sats);              % Reflection status (1 if captured)
sp_stat_arg = zeros(size(sp_stat));           % Reflection over Argentina status (1 if captured)
sp_inc_ang_aux = zeros(N,num_gps_sats);       % Captured reflection incidence angle
sp_inc_ang = cell(num_gps_sats,1);            % Captured reflection incidence angle (cell)
el_sp = zeros(N,num_gps_sats);                % SP elevation angle - RTS
az_sp = zeros(N,num_gps_sats);                % SP azimuth angle - RTS
G_ant_aux = zeros(N,num_gps_sats);            % Antenna gain corresponding to each reflection captured
G_ant = cell(num_gps_sats,1);                 % Antenna gain corresponding to each reflection captured (cell)
sp_path_loss_aux = zeros(N,num_gps_sats);     % Path loss corresponding to each reflection captured
sp_path_loss = cell(num_gps_sats,1);          % Path loss corresponding to each reflection captured (cell)
sp_pres = max(sp_stat,[],2);                  % Reflections present: 1 if at least 1 reflection is captured
sp_time = zeros(1,length(sp_pres)-reg_len_samples);   % Time with reflection presents