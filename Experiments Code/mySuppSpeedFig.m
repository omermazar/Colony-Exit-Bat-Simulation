function [stats] = mySuppSpeedFig(TableCaveExit)

groupBy = {'NumberOfBats', 'VelocityFree' };
lineBy = 'NumberOfBats';
checkParams = {'obsTotalInterferenceRatio', 'dist2Consps_Mean'};
uLines = [5, 10, 40, 100];
cLines = lines(numel(uLines));

%% Grouping
[g, TID]= findgroups(TableCaveExit(:, groupBy));
TableCaveExit.group = g;
analysisTable = TID;
analysisTable.numel  = splitapply(@numel, TableCaveExit.NumberOfBats, g);


%% mean STD GLM
formParams = ' ~ 1 + NumberOfBats + VelocityFree';

for iParam = 1:numel(checkParams)
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
        end % switch
        stats.mdl.(currParam) = fitglm(TableCaveExit, form, 'Distribution', dist);
    end % if ismember(currParam, TableCaveExit.Properties.VariableNames)

end % for iParam


%% Plot
figure; tiledlayout("vertical"); 

for kLine = 1:numel(uLines)
    ax{kLine} = subplot(4,1,kLine); hold on; grid on

    currValue = uLines(kLine);

    ix = true(size(analysisTable,1),1);
    tt = analysisTable(ix,:);

    ixg = ismember(tt.(lineBy), currValue );
    if isnumeric(currValue)
        txtVlue = num2str(currValue);
    else
        txtVlue = char(currValue);
    end
    %%% errorbar of jamming prob
    yyaxis left %
    errorbar(tt.VelocityFree(ixg), tt.obsTotalInterferenceRatio_mean(ixg), ...
        tt.obsTotalInterferenceRatio_std(ixg)  ./ sqrt(tt.numel(ixg)-1), ...
        'o-', 'LineWidth',1.5, DisplayName= [lineBy, txtVlue]) % 'color', cLines(kLine,:)
    
    %%% errorbar of distance
    yyaxis right %
    errorbar(tt.VelocityFree(ixg), tt.dist2Consps_Mean_mean(ixg), ...
        tt.dist2Consps_Mean_std(ixg)  ./ sqrt(tt.numel(ixg)-1), ...
        'o--', 'LineWidth',1.5, DisplayName= [lineBy, txtVlue]) % 'color', cLines(kLine,:)
    % title([num2str(currValue),' bats'], 'FontSize', 9)
    if kLine < numel(uLines)
        xticklabels('')
    end % if kLine
    %%%% Size
    ax{kLine}.Position(2) = 0.08 + (numel(uLines)-kLine)*0.23;
    ax{kLine}.Position(4) = 0.2;
end % for kLine = 1:numel(uLines)

% xticks([3:1:10])
% xticklabels([3:1:10])
xlabel('Flight Speed (m/s)')

stats.analysisTable = analysisTable;



