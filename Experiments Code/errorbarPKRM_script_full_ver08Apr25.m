% plot the rror bar of the PK_RM with and without Masking


%% input
clear TableCaveExit1; % avoid Confuison 
% load pk table drag from:
% C:\Users\YossiYNB3\Documents\University\Research\experiments\caveExit\Figures\forFinal\PK_nBats
tableExperimentPK = TableCaveExit1; 
tableExperimentPK.BatSpecies = repmat(categorical("PK"), height(tableExperimentPK), 1);

clear TableCaveExit1; % avoid Confuison 
% load RM table from:
 % C:\Users\YossiYNB3\Documents\University\Research\experiments\caveExit\Figures\forFinal\Rhino_nBats
tableExperimentRM = TableCaveExit1;
tableExperimentRM.BatSpecies = repmat(categorical("RM"), height(tableExperimentRM), 1);

%% Concatenate the tables
table_PK_RM  = [tableExperimentPK; tableExperimentRM];
% add the differnce criterion
table_PK_RM.PK_plus_Mask = zeros( height(table_PK_RM), 1);
table_PK_RM.PK_plus_Mask((table_PK_RM.BatSpecies == "PK") &  table_PK_RM.MaskingByConsps == 0) = 0;
table_PK_RM.PK_plus_Mask((table_PK_RM.BatSpecies == "PK") &  table_PK_RM.MaskingByConsps == 1) = 1;
table_PK_RM.PK_plus_Mask((table_PK_RM.BatSpecies == "RM") &  table_PK_RM.MaskingByConsps == 0) = 2;
table_PK_RM.PK_plus_Mask((table_PK_RM.BatSpecies == "RM") &  table_PK_RM.MaskingByConsps == 1) = 3;

% remove RM_NoMask from Table (plot only for pk)

% table_PK_RM(table_PK_RM.PK_plus_Mask == 2, :) = [];

%% Stats and plots
groupBy =  {'NumberOfBats', 'MaskingByConsps', 'BatSpecies', 'PK_plus_Mask'};
lineBy  = 'PK_plus_Mask';
% checkParams = { 'ExitSuccess', 'ExitTimesSec', 'CrushesObsPerSec', 'obsTotalInterferenceRatio', 'obsDetectProb_1m'}; %, original
%     % 'NumberOfPulses', 'TotalCaveExitPulses', 'TotalSearchPulses', 'TotalFollowWallPulses', 'TotalBuzzPulses'};
    
checkParams = { 'ExitSuccess', 'ExitTimesSec', 'CrushesObsPerSec', 'obsTotalInterferenceRatio', 'obsDetectProb_1m', ...
    'NumberOfPulses', 'TotalCaveExitPulses', 'TotalSearchPulses', 'TotalFollowWallPulses', 'TotalBuzzPulses', ...
    'TotalObstacleManuverPulses', 'TotalAvoidBatManPulses'}; % new

[statsPKRM, table_full] = analyzeCaveExitExpNew(table_PK_RM, groupBy, lineBy, 'none', checkParams);

%%
% stats_PK = load('C:\Users\YossiYNB3\Documents\University\Research\experiments\caveExit\Figures\forFinal\PK_nBats\Stats.mat');
% stats_RM = load('C:\Users\YossiYNB3\Documents\University\Research\experiments\caveExit\Figures\forFinal\Rhino_nBats\Stats.mat');
% stats_PK_RM =  load('C:\Users\YossiYNB3\Documents\University\Research\experiments\caveExit\Figures\forFinal\PK_RM_Comp\Stats.mat');
% 
% %% Concatenate the tables
% % pk
% table_PK = stats_PK.analysisTable;
% % add the batSpeciest
% table_PK.BatSpecies = repmat(categorical("PK"), height(table_PK),1);
% table_PK = movevars(table_PK, "BatSpecies", "After", "NumberOfBats");

% % RM
% table_RM = stats_RM.analysisTable;
% % add the batSpeciest
% table_RM.BatSpecies = repmat(categorical("RM"), height(table_RM),1);
% table_RM = movevars(table_RM, "BatSpecies", "After", "NumberOfBats");
% table_RM.ObsMemorySize = [];
% 
% table_PK_RM = [table_PK; table_RM];

%% the plot
tt = statsPKRM.analysisTable;
params = ["ExitSuccess", "ExitTimesSec", "CrushesObsPerSec", "obsTotalInterferenceRatio", "obsDetectProb_1m"]; % original
% params = ["NumberOfPulses"]; % new
% params = string(checkParams);
% iParam = "ExitSuccess"; 
for kFig = 1:numel(params)
    iParam = params(kFig);
    figure; hold on; grid on
    meanName   = strcat(iParam,'_mean');
    stdName    = strcat(iParam,'_std');
    uVals = unique(table_PK_RM.PK_plus_Mask)';
    
    line_style = {"--", "-", "--", "-"};
    color_line = {'k', 'k', [ 0.7 0.7 0.7], [0.7 0.7 0.7]};
    for kVal = uVals
        ixg = tt.PK_plus_Mask == kVal;
        if ismember(kVal, [0,1]) % PK
            color_line = 'k';
        elseif ismember(kVal, [2,3]) % RM
            color_line = [0.6, 0.6, 0.6];
        end
        if ismember(kVal, [0,2]) % No Masking
            line_style = '--';
        elseif ismember(kVal, [1,3]) % Masking
            line_style = '-';
        end
        iLine= find(uVals == kVal);

        errorbar(tt.NumberOfBats(ixg), tt.(meanName)(ixg), tt.(stdName)(ixg)  ./ sqrt(tt.numel(ixg)-1), ...
            'o', 'LineWidth',1.5, 'lineStyle',  line_style, 'Color', color_line, 'DisplayName', num2str(kVal))
        xlabel('Bat number')
        xticks([0:20:max(xlim)])
    end % for kVal
    title(iParam)
end % for kFig

  


