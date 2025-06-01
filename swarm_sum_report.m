function [swarm_strcut_full] = swarm_sum_report(BatDATA,Percentile, pick_mid_bats_flag)

% summary of the trial
if nargin < 3
    pick_mid_bats_flag = false; % pick 25-percentile bats to 75% only and divide the x-range by two
end % if nargin <3

if nargin < 2
   Percentile = [50 80 90]; 
end % naragi

%% summary
swarm_strcut_full.num_of_bats = BatDATA.AllParams.SimParams.TotalBatsNumber;
swarm_strcut_full.run_time = BatDATA.AllParams.SimParams.SimulationTime;
swarm_strcut_full.yDist_std = BatDATA.AllParams.SimParams.swarm_initial_ybat_std;
swarm_strcut_full.xRange = BatDATA.AllParams.SimParams.swarm_initial_xbat_max;
swarm_strcut_full.start_angle_std = rad2deg(BatDATA.AllParams.SimParams.swarm_initial_teta_std);
swarm_strcut_full.receiver_type = BatDATA.AllParams.BatSonarParams.ReceiverType;
swarm_strcut_full.flight_type = BatDATA.AllParams.BatFlightParams.FlightTypeFlag;

if pick_mid_bats_flag
   swarm_strcut_full.num_of_bats = round(swarm_strcut_full.num_of_bats / 2);
   swarm_strcut_full.xRange = swarm_strcut_full.xRange / 2;
end

%% Distance from cave
% load('D:\Dropbox\University\מחקר\experiments\Swarm\Nvigation DATA\DistFromCave_distFromMidSwarm.mat');
% distance from cave = xBin_distRangeCave;
dist_cave_in = [50,150,250,350,450,550,650,750,850,950,1050,1150,1250,1350,1450,1550,1650,1750,1850,1950];
dist_cave_in = [0, dist_cave_in, 4400];
% mean of distnce to mid-swarm = dataDistMidSwarm_all_smooth
mean_yin = [9.68,20.91,28.42,31.81,35.21,37.59,37.13,41.32,50.25,58.77,63.77,69.64,76.15,78.74,77.47,81.34,84.32,85.66,90.45,88.32];
std_y = sqrt(pi/2)*mean_yin;
std_y = [2.5, std_y, 250];
% Interpolation of the distance
swarm_strcut_full.xDist_2Cave = interp1(std_y, dist_cave_in, swarm_strcut_full.yDist_std)';

% % Distance from cace approximation from the fiעures of swarm 
% range_from_cave_in = [ 1 2 250 750 1250 1750]; % approximation from the figure
% dist_drom_mid_in = [5 10 75 150 200 250];
% Interpolation of the distance
% swarm_strcut_full.xDist_2Cave = interp1(dist_drom_mid_in, range_from_cave_in, swarm_strcut_full.yDist_std);

%% summarize
swarm_strcut_full.AllParams = BatDATA.AllParams;

swarm_strcut_full.interference_analysis = BatDATA.SwarmSummary;

swarm_strcut_full.percentile_bats = find_prctile(BatDATA, Percentile);

swarm_strcut_full.all_bats.BAT(swarm_strcut_full.num_of_bats) = ...
    struct('BatNum',[], ....
    'start_pos', [], ...
    'start_angle', [], ...
    'xBat', [], ...
    'yBat', [], ...
    'own_calls_lvl', [], ...
    'rx_lvl_from_conspecific', [], ...
    'masking_lvl', [], ...
    'consps_in_each_calls', [] ...
    );
%%%% Pick the mid-bats
if ~pick_mid_bats_flag
    idx_pick_bats = 1:1:swarm_strcut_full.num_of_bats;
else % if ~pick_mid_bats_flag
    low_prc = 25; high_prc = 75;
    all_x0 = [BatDATA.BAT.BatX0];
    [~,q] = sort(all_x0);
    n = numel(all_x0);
    first_bat_ix = round(n*low_prc/100);
    last_bat_ix = round(n*high_prc/100)-1;
    idx_pick_bats = q(first_bat_ix:1:last_bat_ix);
    % if there are less bats than expected add to the end of the swarm
    if numel(idx_pick_bats) < swarm_strcut_full.num_of_bats
        ndiff = swarm_strcut_full.num_of_bats - numel(idx_pick_bats);
        idx_pick_bats = [idx_pick_bats, q(last_bat_ix:last_bat_ix+ndiff)];
    end % if numel
%     ix_mid_swarm = all_x0 <= prctile(all_x0, high_prc) & all_x0 >=  prctile(all_x0, low_prc)
end % if ~pick_mid_bats_flag
BAT = BatDATA.BAT;
start_y = [BAT.BatY0];

for kBat = 1:numel(idx_pick_bats)
    orig_bat = idx_pick_bats(kBat);
    swarm_strcut_full.all_bats.BAT(kBat).BatNum = orig_bat;
    swarm_strcut_full.all_bats.BAT(kBat).start_pos = [BAT(orig_bat).BatX0, BAT(orig_bat).BatY0];
    swarm_strcut_full.all_bats.BAT(kBat).start_angle = BAT(orig_bat).Teta(1);
    swarm_strcut_full.all_bats.BAT(kBat).start_percentile = sum(start_y <= BAT(orig_bat).BatY0 ) / numel(start_y);
    swarm_strcut_full.all_bats.BAT(kBat).xBat = BAT(orig_bat).xBatPos;
    swarm_strcut_full.all_bats.BAT(kBat).yBat = BAT(orig_bat).yBatPos;
    swarm_strcut_full.all_bats.BAT(kBat).own_calls_lvl = BAT(orig_bat).BatSonarEchosMat(1,:);
    swarm_strcut_full.all_bats.BAT(kBat).rx_lvl_from_conspecific = BAT(orig_bat).Consps_EchosVec;
    swarm_strcut_full.all_bats.BAT(kBat).masking_lvl = BAT(orig_bat).AllInterPulses;
    swarm_strcut_full.all_bats.BAT(kBat).consps_in_each_calls =  BAT(orig_bat).Consps_FindsStruct;
end %for kBat

    
    
