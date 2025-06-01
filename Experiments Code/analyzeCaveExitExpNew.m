function [stats, TableCaveExit] = analyzeCaveExitExpNew(TableCaveExit, groupBy, lineBy, subplotBy, checkParams, saveResFlg, savePath)

if nargin < 2
%     groupBy = {'NumberOfBats', 'MaskingByConsps', 'ClutterGainDB' ,
%     'ObsMemorySize', 'PulsePower', 'VelocityManuever'}; % 'ConfusionObs'
    groupBy = {'NumberOfBats', 'MaskingByConsps', 'BatSpecies' }; 
end

if nargin < 3
   lineBy = 'MaskingByConsps'; % 'ObsMemorySize'; % 'ClutterGainDB'; 'PulsePower; ;'MaskingByConsps', 'ConfusionObs'
end

if nargin < 4
    subplotBy = 'none'; %  'none' or any fields 
end

if nargin < 5
    noCheck = true;
elseif isempty(checkParams)
    noCheck = true;
else
    noCheck = false;
end
if noCheck
    % checkParams = {'ExitSuccess', 'ExitTimesSec', 'CrushesObsPerSec', ...
    %     'obsTotalInterferenceRatio', 'conspsTotalInterferenceRatio', ...
    %     'obsDetectProb_1m', 'conspsDetectProb_1m', 'CrushesConspsPerSec', 'CrushesConspsTotal', 'CrushesObsTotal'};
    checkParams = {'ExitSuccess', 'ExitTimesSec', 'CrushesObsPerSec', 'obsTotalInterferenceRatio'};
end

if nargin < 6
    saveResFlg = 0;
end

%% Filter Out DATA
% remove rows of detect only
% ixDetOnly = TableCaveExit.MaskingByConsps == 0 & TableCaveExit.DetectConsps == 1;
% TableCaveExit = TableCaveExit(~ixDetOnly,:);

%% Calcualte rate of parameters and add to table 

if ~ismember('CrushesObsPerSec', TableCaveExit.Properties.VariableNames) && ~any(ismember(TableCaveExit.Properties.VariableNames, 'TotalFlightTime')) 
    TableCaveExit.TotalFlightTime = TableCaveExit.ExitTimesSec;
    TableCaveExit.TotalFlightTime(isnan(TableCaveExit.ExitTimesSec)) = 15;
    TableCaveExit.CrushesObsPerSec = TableCaveExit.CrushesObsTotal ./ TableCaveExit.TotalFlightTime;
    TableCaveExit.CrushesConspsPerSec = TableCaveExit.CrushesConspsTotal ./ TableCaveExit.TotalFlightTime;
end % if ~any

%% Add bat species collumn 
if ~any(ismember(TableCaveExit.Properties.VariableNames, 'BatSpecies'))
    TableCaveExit.BatSpecies = repmat(categorical("PK"), size(TableCaveExit,1), 1);
end % if ~any
%% grouping
% group by bats number
% groupBy = 'NumberOfBats';
% [gr, nBats] = findgroups(TableCaveExit.(groupBy));

plotFlag = true;
% mean and std of the measured params
% checkParams = {'ExitSuccess', 'ExitTimesSec', 'CrushesObsTotal', 'CrushesConspsTotal', 'TotalFilghtDistance', ...
%     'TotalObstacleManuverPulses' 'TotalObstacleManuverSequences', ...
%     'TotalAvoidBatManPulses', 'TotalAvoidBatManSequences', 'TotalCaveExitPulses', 'TotalCaveExitSequences', ...
%     'TotalFollowWallSequences', 'TotalFollowWallPulses', ...
%     'obsTotalInterferenceRatio', 'conspsTotalInterferenceRatio'};
% checkParams = {'ExitSuccess', 'ExitTimesSec', 'TotalFilghtDistance', 'CrushesObsTotal', 'CrushesConspsTotal', ...
%     'dist2Consps_Med', 'dist2Consps_Min', 'obsTotalInterferenceRatio', 'conspsTotalInterferenceRatio''obsDetectProb_1m', 'obsDetectProb_2m', 'obsDetectProb_2m'};
% checkParams = {'TotalFilghtDistance', 'dist2Consps_Med', 'TotalFilghtDistance', 'dist2Consps_Min'};

%%%%


%%%%% for special Tests
% checkParams = {'ExitSuccess', 'ExitTimesSec', 'CrushesObsPerSec', ...
%      'obsTotalInterferenceRatio', 'conspsTotalInterferenceRatio', ...
%      'obsDetectProb_1m', 'conspsDetectProb_1m'};
% checkParams = {'dist2Consps_Mean', 'dist2Consps_Med', 'dist2Consps_Min'};
% checkParams = {'TotalObstacleManuverPulses', 'TotalAvoidBatManPulses', 'TotalCaveExitPulses', 'TotalFollowWallPulses', ...
    % 'conspsTotalInterferenceRatio'};
% { 'obsDetectProb_1m', 'conspsDetectProb_1m', 'obsDetectProb_2m', 'conspsDetectProb_2m', 'obsDetectProb_3m', 'conspsDetectProb_3m',};
% checkParams = {'ExitSuccess', 'ExitTimesSec', 'CrushesObsPerSec', ...
%      'obsTotalInterferenceRatio', 'conspsTotalInterferenceRatio', ...
%      'obsDetectProb_1m', 'conspsDetectProb_1m'};
%%%

% testParam = 'ClutterGainDB'; %'ObsMemorySize'; % 'ClutterGainDB'; 'PulsePower; 'none' ;'MaskingByConsps''; 'VelocityManuever'
%%%%

%% new 2024- Add Cluser and Obs 
if strcmp(lineBy, 'ExpTypeConfCluster') && ~ismember('ExpTypeConfCluster', TableCaveExit.Properties.VariableNames)
    TableCaveExit.ExpTypeConfCluster = TableCaveExit.ConfusionObs + TableCaveExit.ObsClusterFlag; 
end %
if ~ismember('ExpTypeConfCluster', groupBy) && ismember('ExpTypeConfCluster', TableCaveExit.Properties.VariableNames)
    groupBy{end+1} = 'ExpTypeConfCluster';
end
%%
% summary Data
[g, TID]= findgroups(TableCaveExit(:, groupBy));
TableCaveExit.group = g;

% [g, batsNum, masking, detConsps]  = findgroups(TableCaveExit.NumberOfBats, TableCaveExit.MaskingByConsps, TableCaveExit.DetectConsps);

stats.summaryData.totalExperiments = size(TableCaveExit,1);
stats.summaryData.numOfOBats       = unique(TableCaveExit.NumberOfBats)';
stats.summaryData.MaskingByConsps  = unique(TableCaveExit.MaskingByConsps)';
if ismember('DetectConsps', TableCaveExit.Properties.VariableNames)
    stats.summaryData.DetectConsps     = unique(TableCaveExit.DetectConsps)';
else
    stats.summaryData.DetectConsps = true(size(stats.summaryData.MaskingByConsps));
end
analysisTable = TID;
analysisTable.numel  = splitapply(@numel, TableCaveExit.NumberOfBats, g);

if strcmp(lineBy, 'NumberOfBats')
    uLines = [1, 5, 10, 40, 100]; % [2, 5, 10, 40]; % plot only 4 lines
    % uLines = unique(analysisTable.(lineBy)); % plot all bat-densities
else
    uLines = unique(analysisTable.(lineBy));
end % if strcmp
% uLines = unique(analysisTable.(lineBy));
cLines = lines(numel(uLines));
% cLines = [0 0 0; 0.6 0.6 0.6];

switch subplotBy
    case 'none'
        nSplits  = 1;
        titleTxt =  '';
        uSplits = '';
        formParams = [' ~ 1 + NumberOfBats + ', groupBy{end} ] ;   % ' ~ 1 + NumberOfBats + MaskingByConsps ' 
        % formParams = [' ~ 1 + ', groupBy{end} ] ;
    otherwise
        uSplits    = unique(analysisTable.(subplotBy));
        nSplits    = numel(unique(uSplits));
        titleTxt   =  [subplotBy, '='];
        formParams = [' ~ 1 + NumberOfBats + MaskingByConsps + ', groupBy{end} ];
end

for iParam = 1:numel(checkParams)
    %% mean and STD
    currParam = checkParams{iParam};
    if ismember(currParam, TableCaveExit.Properties.VariableNames)
        meanName   = [currParam,'_mean'];
        stdName    = [currParam,'_std'];
        analysisTable.(meanName)  = splitapply(@nanmean, TableCaveExit.(currParam), g) ;
        analysisTable.(stdName)   = splitapply(@nanstd , TableCaveExit.(currParam), g) ;

        %% glm
        % linear model
        form = [currParam, formParams];
        switch currParam
            case {'ExitSuccess', 'obsTotalInterferenceRatio', 'conspsTotalInterferenceRatio', ...
                'obsDetectProb_1m', 'conspsDetectProb_1m'} 
                dist = 'binomial';
            case {'CrushesObsPerSec', 'CrushesConspsPerSec', 'CrushesObsTotal', 'CrushesConspsTotal' }
                dist = 'poisson';
            case 'ExitTimesSec'
                dist = 'normal';
            otherwise
                dist = 'normal';
        end % if ismember
       
        try
            stats.mdl.(currParam) = fitglm(TableCaveExit, form, 'Distribution', dist);
            form2 = [form, ' + ', groupBy{end}, '^2'];
            stats.mdl.([currParam,'Square']) = fitglm(TableCaveExit, form2, 'Distribution', dist);
        catch
            warning(strcat(form, " ... was not calculated!!!"))
        end % try

        %% Plot
        if plotFlag
            fig = figure; hold on;
            for kSplit = 1:nSplits
                ax{kSplit} = subplot(1, nSplits, kSplit); hold on
                % The groups to Plot
                %         ixMasking    = analysisTable.masking == 1 & analysisTable.detConsps == 1;
                % %         ixDetectOnly = analysisTable.masking == 0 & analysisTable.detConsps == 1;
                %         ixNoMasking  = analysisTable.masking == 0 & analysisTable.detConsps == 0;
                % filter the table vy the splits
                if ~strcmp(subplotBy, 'none')
                    ix = ismember(analysisTable.(subplotBy), uSplits(kSplit));
                else
                    ix = true(size(analysisTable,1),1);
                end
                tt = analysisTable(ix,:);
                for kLine = 1:numel(uLines)
                    currValue = uLines(kLine);
                    ixg = ismember(tt.(lineBy), currValue );
                    if isnumeric(currValue)
                        txtVlue = num2str(currValue);
                    elseif islogical(currValue)
                        txtVlue = num2str(currValue);
                    else
                        txtVlue = char(currValue);
                    end
                    if ~strcmp(lineBy, 'NumberOfBats')
                        errorbar(ax{kSplit}, tt.NumberOfBats(ixg), tt.(meanName)(ixg), ...
                            tt.(stdName)(ixg)  ./ sqrt(tt.numel(ixg)-1), ...
                            'o-', 'color', cLines(kLine,:), 'LineWidth',1.5, DisplayName= [lineBy, txtVlue])
                         xlabel('Bat number')
                         xticks([0:20:max(xlim)])
                    else
                        %%% errorbar as function of measured params
                        errorbar(ax{kSplit}, tt.(groupBy{end})(ixg), tt.(meanName)(ixg), ...
                            tt.(stdName)(ixg)  ./ sqrt(tt.numel(ixg)-1), ...
                            'o-', 'color', cLines(kLine,:), 'LineWidth',1.5, DisplayName= [lineBy, txtVlue])
                         xlabel(groupBy{end})
%                          xticks(0:1:max(xlim))
                    end % if ~strcmp(lineBy, 'NumberOfBats')
                end % for kLine = 1:numel(uLines)
                
               
                %             title([ splitBy, '=', num2str(uSplits(kSplit))], 'Interpreter','none')
                if ~strcmp(subplotBy, 'none')
                    title([ titleTxt, num2str(uSplits(kSplit))], 'Interpreter','none')
                end
                grid
                if kSplit ==1
                    L = legend;
                    L.Box = 'off';
                    L.Location = 'best';
                end % if kSplit
                %         fig.Position= [100 100 280 250];
            end % for kax

            %%%% Final Position and title
            linkaxes
            if nSplits >1
                fig.Position = [100 300 450*nSplits 420];
            else
                fig.Position = [100 300 250*nSplits 420];
            end
            sgtitle(checkParams{iParam},'Interpreter', 'none')%            
          
            %% Saving the figures
            if saveResFlg
                % final for save
                fileName = ['CaveExit_CLutter_', currParam];
                savefig(fig, fullfile(savePath, fileName))
                
                xlabel('')

                set(fig,'Units','inches')
                fig.Position([3,4]) = [2.5, 4]/1.1; % [1.5, 1.9]/1.1;% for small panels 
                sgtitle('')
                ylabel('')
                legend('off')
                xlabel('')               
                print(fig, fullfile(savePath, fileName), '-djpeg', '-r100')
            end % if saveFig

        end % if plotFlag
    end % if ismember(currParam, TableCaveExit.Properties.VariableNames)
end % for iParam
            

stats.analysisTable = analysisTable;

if saveResFlg
    save(fullfile(savePath, 'Stats.mat'), 'stats')
    writetable(TableCaveExit, fullfile(savePath, 'TableCaveExit_New_intWin_.txt'))
end

