function [stats] = analyzeCaveExitExp(TableCaveExit, groupBy, saveResFlg, savePath)

if nargin < 2
%     groupBy = {'NumberOfBats', 'MaskingByConsps', 'ClutterGainDB' , 'ObsMemorySize', 'PulsePower'};
    groupBy = {'NumberOfBats', 'MaskingByConsps', 'ClutterGainDB' }; 
end

if nargin < 4
    saveResFlg = 0;
end

%% Filter Out DATA
% remove rows of detect only
ixDetOnly = TableCaveExit.MaskingByConsps == 0 & TableCaveExit.DetectConsps == 1;
TableCaveExit = TableCaveExit(~ixDetOnly,:);

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
checkParams = {'ExitSuccess', 'ExitTimesSec', 'CrushesObsTotal', 'CrushesConspsTotal', ...
     'obsTotalInterferenceRatio', 'conspsTotalInterferenceRatio', ...
     'obsDetectProb_1m', 'obsDetectProb_2m', 'obsDetectProb_3m', ...
     'conspsDetectProb_1m', 'conspsDetectProb_2m', 'conspsDetectProb_3m'};
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

gr1 = 'none'; % 'ObsMemorySize'; % 'ClutterGainDB'; 'PulsePower; 'none' 'MaskingByConsps'
gr2 = 'MaskingByConsps';
axesBy = 'gr1'; % 'gr1' | 'gr2'

switch axesBy
    case 'gr2'
        lineBy    = gr1;
        subplotBy = gr2;
        fN = 'byParam';
        lTxt = gr1(1:7);
    case 'gr1'
        lineBy    = gr2;
        subplotBy = gr1;
        fN = 'byMasking';
        lTxt = 'masking';
end % switch axesBy
uLines = unique(analysisTable.(lineBy));
cLines = lines(numel(uLines));

switch subplotBy
    case 'none'
        nSplits  = 1;
        titleTxt =  '';
        uSplits = '';
        formParams = ' ~ 1 + NumberOfBats + MaskingByConsps ' ;
    otherwise
        uSplits    = unique(analysisTable.(subplotBy));
        nSplits    = numel(unique(uSplits));
        titleTxt   =  [subplotBy, '='];
        formParams = [' ~ 1 + NumberOfBats + MaskingByConsps + ', gr1 ];
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

                    errorbar(ax{kSplit}, tt.NumberOfBats(ixg), tt.(meanName)(ixg), ...
                        tt.(meanName)(ixg)  ./ sqrt(tt.numel(ixg)-1), ...
                        'o-', 'color', cLines(kLine,:), 'LineWidth',1.5, DisplayName= [lTxt, num2str(currValue)])
                end % for kLine = 1:numel(uLines)
                % NoMasking
                %             errorbar(batsNum(ixNoMasking), analysisTable.(meanName)(ixNoMasking), ...
                %                 analysisTable.(meanName)(ixNoMasking)  ./ sqrt(analysisTable.numel(ixNoMasking)-1), ...
                %                 'o-k', 'LineWidth',1.5, DisplayName='None')

                % DetectOnly
                %          errorbar(batsNum(ixDetectOnly), analysisTable.(meanName)(ixDetectOnly), ...
                %              analysisTable.(meanName)(ixDetectOnly)  ./ sqrt(analysisTable.numel(ixDetectOnly)-1), ...
                %             'o-b', 'LineWidth',1.5, DisplayName='DetectOnly')

                xlabel('Bats number')
                %             title([ splitBy, '=', num2str(uSplits(kSplit))], 'Interpreter','none')
                if ~strcmp(subplotBy, 'none')
                    title([ titleTxt, num2str(uSplits(kSplit))], 'Interpreter','none')
                end
                grid
                if kSplit ==1
                    L = legend;
                    L.Box = 'off';
                end % if kSplit
                %         fig.Position= [100 100 280 250];
            end % for kax

            %%%% Final Position and title
            linkaxes
            fig.Position = [100 300 250*nSplits 420];
            sgtitle(checkParams{iParam},'Interpreter', 'none')

            %% Saving the figures
            if saveResFlg
                fileName = ['caveExit_', fN, '_', currParam];
                savefig(fig, fullfile(savePath, fileName))
                print(fig, fullfile(savePath, fileName), '-djpeg', '-r300')
            end % if saveFig

        end % if plotFlag
    end % if ismember(currParam, TableCaveExit.Properties.VariableNames)
end % for iParam
            

stats.analysisTable = analysisTable;

if saveResFlg
    save(fullfile(savePath, 'Stats.mat'), 'stats')
    writetable(TableCaveExit, fullfile(savePath, 'TableCaveExit_1.txt'))
end

