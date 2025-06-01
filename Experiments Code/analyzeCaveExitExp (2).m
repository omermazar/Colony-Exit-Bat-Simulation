function [stats] = analyzeCaveExitExp(TableCaveExit, groupBy, lineBy, subplotBy, saveResFlg, savePath)

if nargin < 2
%     groupBy = {'NumberOfBats', 'MaskingByConsps', 'ClutterGainDB' , 'ObsMemorySize', 'PulsePower', 'VelocityManuever'};
    groupBy = {'NumberOfBats', 'MaskingByConsps', 'ClutterGainDB' }; 
end

if nargin < 3
   lineBy = 'MaskingByConsps'; % 'ObsMemorySize'; % 'ClutterGainDB'; 'PulsePower; ;'MaskingByConsps' 
end

if nargin < 4
    subplotBy = 'none'; %  'none' or any fields 
end

if nargin < 5
    saveResFlg = 0;
end

%% Filter Out DATA
% remove rows of detect only
% ixDetOnly = TableCaveExit.MaskingByConsps == 0 & TableCaveExit.DetectConsps == 1;
% TableCaveExit = TableCaveExit(~ixDetOnly,:);

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
checkParams = {'ExitSuccess', 'ExitTimesSec', 'CrushesObsTotal', ...
     'obsTotalInterferenceRatio', 'obsManueverFalseAlarm_ratio', ...
     'obsDetectProb_1m', 'TotalFilghtDistance'};
%%%%% for special Tests
%  checkParams = {'conspsDetectProb_1m', 'conspsTotalInterferenceRatio', 'CrushesConspsTotal' }
% checkParams = {'ExitSuccess', 'ExitTimesSec', 'CrushesObsTotal', ...
%      'obsTotalInterferenceRatio',...
%      'obsDetectProb_1m', 'TotalFilghtDistance', ...
%      'CrushesConspsTotal', 'conspsTotalInterferenceRatio' , 'conspsDetectProb_1m', ...
%      'TotalAvoidBatManSequences', 'TotalAvoidBatManPulses',};
%  checkParams = {'conspsDetectProb_1m', 'CrushesConspsTotal', 'TotalSearchPulses', 'TotalAvoidBatManPulses', ...
%      'TotalCaveExitPulses', 'conspsTotalInterferenceRatio', 'TotalFollowWallPulses', 'TotalObstacleManuverPulses', ...
%      'TotalCaveExitSequences', 'TotalFollowWallSequences'  'TotalObstacleManuverSequences'};
% %  checkParams = {'obsManueverFalseAlarm_ratio'}; % {'TotalAvoidBatManPulses'}; %{'CrushesConspsTotal'};
%{'TotalObstacleManuverPulses'}; % { 'obsDetectProb_1m', 'conspsDetectProb_1m', 'obsDetectProb_2m', 'conspsDetectProb_2m', 'obsDetectProb_3m', 'conspsDetectProb_3m',};
%%%

% testParam = 'VelocityManuever'; %'ObsMemorySize'; % 'ClutterGainDB'; 'PulsePower; 'none' ;'MaskingByConsps''; 'VelocityManuever'
%%%%

% summary Data
[g, TID]= findgroups(TableCaveExit(:, groupBy));
TableCaveExit.group = g;

% [g, batsNum, masking, detConsps]  = findgroups(TableCaveExit.NumberOfBats, TableCaveExit.MaskingByConsps, TableCaveExit.DetectConsps);

stats.summaryData.totalExperiments = size(TableCaveExit,1);
stats.summaryData.numOfOBats       = unique(TableCaveExit.NumberOfBats)';
stats.summaryData.MaskingByConsps  = unique(TableCaveExit.MaskingByConsps)';
stats.summaryData.DetectConsps     = unique(TableCaveExit.DetectConsps)';
analysisTable = TID;
analysisTable.numel  = splitapply(@numel, TableCaveExit.NumberOfBats, g);

if strcmp(lineBy, 'NumberOfBats')
    uLines = [1,5,40,100];
else % if strcmp(lineBy, 'NumberOfBats')
    uLines = unique(analysisTable.(lineBy));
end % if strcmp(lineBy, 'NumberOfBats')
cLines = lines(numel(uLines));
% cLines = [0 0 0; 0.6 0.6 0.6];

switch subplotBy
    case 'none'
        nSplits  = 1;
        titleTxt =  '';
        uSplits = '';
        formParams = [' ~ 1 + NumberOfBats + ', groupBy{end} ] ;   % ' ~ 1 + NumberOfBats + MaskingByConsps ' 
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
        %     if islogical(TableCaveExit.(currParam))
        %         TableCaveExit.(currParam) = double(TableCaveExit.(currParam ));
        %     end
        stats.mdl.(currParam) = fitglm(TableCaveExit, form);

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
                    else
                        txtVlue = char(currValue);
                    end
                    if ~strcmp(lineBy, 'NumberOfBats')
                        errorbar(ax{kSplit}, tt.NumberOfBats(ixg), tt.(meanName)(ixg), ...
                            tt.(stdName)(ixg)  ./ sqrt(tt.numel(ixg)-1), ...
                            'o-', 'color', cLines(kLine,:), 'LineWidth',1.5, 'DisplayName', [lineBy, txtVlue])
                         xlabel('Bat number')
                         xticks([0:10:max(xlim)])
                    else
                        %%% errorbar as function of measured params
                        errorbar(ax{kSplit}, tt.(groupBy{end})(ixg), tt.(meanName)(ixg), ...
                            tt.(stdName)(ixg)  ./ sqrt(tt.numel(ixg)-1), ...
                            'o-', 'color', cLines(kLine,:), 'LineWidth',1.5, 'DisplayName', [lineBy, txtVlue])
                         xlabel(groupBy{end})
                         xticks(0:2:max(xlim))
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
                fileName = ['caveExit_Final_MemSize_', currParam];
                savefig(fig, fullfile(savePath, fileName))

                set(fig,'Units','inches')
                fig.Position([3,4]) = [2.5, 4]/1.1; % [1.5, 1.9]/1.1;% for small panels 
                sgtitle('')
                ylabel('')
                L.delete
                xlabel('')
                delete( fig.Children(2))
                print(fig, fullfile(savePath, fileName), '-djpeg', '-r100')
            end % if saveFig

        end % if plotFlag
    end % if ismember(currParam, TableCaveExit.Properties.VariableNames)
end % for iParam
            

stats.analysisTable = analysisTable;
if saveResFlg
    save(fullfile(savePath, 'Stats.mat'), 'stats')
    writetable(TableCaveExit, fullfile(savePath, 'TableCaveExit_MemSize.txt'))
end

