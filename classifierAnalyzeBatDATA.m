function [sumTable, stats] = classifierAnalyzeBatDATA(BatDATA, BatNum)

BAT = BatDATA.BAT(BatNum);
NPulses = BAT.NumOfTimesSonarTransmits;


%% Build Summary Table of all filterBank Detections
pulseResults = struct();
jRow = 1;
sumTable = [];
for kPulse =1:NPulses
    try
        currStruct = BAT.FindsMaskingStruct(kPulse) ;
        % row to column
        if ~isempty(currStruct.ClassifierFullResults)
            fn = fieldnames(currStruct.ClassifierFullResults);
            %     nCols = numel(fn);
            for iField = 1:numel(fn)
                results.(fn{iField}) =  currStruct.ClassifierFullResults.(fn{iField})';
            end % for iField
                
            % Classifier True-Positive
            results.DetectedAndTrue = ismember(results.ntime, currStruct.ClassifierDetecdtedAndTRUE);
            % False Alarms (True Negative)
            results.FalseAlarm = ismember(results.ntime, currStruct.ClassifierFalseAlarmsTimes);
            results.MissByJamm = ismember(results.ntime, currStruct.ClassifierMaskedTimes);
            %     % Classifier Misses (JAmmed)- Not Seen
            %     results.Miss = ismember(results.ntime, currStruct.ClassifierMaskedTimes);
            currTable = struct2table(results);
            
            nCols =  size(currTable,2);

            %% add summary data from the Struct
            %%%% add mis-detections by FilterBank
            if ~isempty(currStruct.CorrelationMissedTargetd)
                % add nans
                nMissDetection = numel(currStruct.CorrelationMissedTargetd);
                nanArray       = [nan(nMissDetection, 3), zeros(nMissDetection, nCols-3)];% nan(nMissDetection, nCols);
                currTable      = [currTable; array2table(nanArray,'variablenames',currTable.Properties.VariableNames)];
                
                % commplete data to new rows
                newRowsIx = (size(currTable,1) - nMissDetection + 1) : size(currTable,1) ;
                
                % Find the right order of the targets
                [~, preDetectedPreys ] = IsPreyDetected('Prey', ...
                    BAT.EchosFromPreyStruct(kPulse), BAT.TransmittedPulsesStruct(kPulse), [], [], BatDATA.AllParams);
                ixDetectedPRey = ismember(preDetectedPreys, currStruct.DetectedPreys);
                ixMissFB = ismember(preDetectedPreys, currStruct.CorrelationMissedTargetd);
            
            else %  if ~isempty(currStruct.CorrelationMissedTargetd)
                newRowsIx = [];
                ixDetectedPRey = true(size(currStruct.DetectedPreys));
                ixMissFB = [];
            end %  if ~isempty(currStruct.CorrelationMissedTargetd)
            nRows = size(currTable,1);
            
            % add more data
            currTable.isTarget = currTable.DetectedAndTrue + currTable.MissByJamm;
            currTable.isTarget(newRowsIx) = ones( numel(newRowsIx), 1); % complete the miss detection
            currTable.MissedDetectedFB = zeros( nRows, 1);
            currTable.MissedDetectedFB(newRowsIx) = ones( numel(newRowsIx), 1); % complete the miss detection
            currTable.pulseNum = repmat(currStruct.PulseNum, nRows,1);
            currTable.callType = repmat(string(BAT.ManueverCmdStruct(kPulse).ManueverStage), nRows,1);
            
            % add Rx-Level by Target
            % Detected Targets are sorted by time of arrival
            currTable.TargetRxLevel = nan(size(currTable.ntime));
            currTable.TargetSIRdB = nan(size(currTable.ntime));
            currTable.TargetRxLevel(ixDetectedPRey) = currStruct.DetectedPreysRxPower(ixDetectedPRey);
            currTable.TargetSIRdB(ixDetectedPRey)   = currStruct.DetecetedPrey2InterferenceRatioDB; % SIR is cmouted only for detecteions
            currTable.TargetRxLevel(newRowsIx)      = currStruct.DetectedPreysRxPower(ixMissFB);
%             currTable.TargetSIRdB(newRowsIx)        = currStruct.DetecetedPrey2InterferenceRatioDB(ixMissFB);

            %%  concatenate
            sumTable = [sumTable; currTable];
        end % if ~isempty(currStruct.ClassifierFullResults)
    catch
        warning(['classifierAnalyzeBatDATA - Pulse: ', num2str(kPulse)]);

    end % try
end % for kPulse

%% Anaze Sum Table
ffn = { 'corrCoeffFR', 'corrCoeffMatrix'};
%%%% {'pks_trghs_rms', 'corrCoeffFR', 'corrCoeffMatrix', 'mseFR', 'mseMatrix'};
nVars = numel(ffn);
nCols = 2;
nRows = ceil(nVars/nCols);

% Stats
stats.TotalFilterBankDetections = sum(~isnan(sumTable.ntime),'omitnan');
stats.TotalTargetsDetections    = sum(sumTable.isTarget,'omitnan');
stats.TotalNonTargets           = sum(~sumTable.isTarget,'omitnan');

stats.nMissesByFB               = sum(sumTable.MissedDetectedFB,'omitnan');   
stats.nClassifiedAsTargets      = sum(sumTable.isClassified,'omitnan');
stats.TruePositive              = sum(sumTable.DetectedAndTrue,'omitnan');
stats.nFalseAlarms              = sum(sumTable.FalseAlarm,'omitnan');
stats.MissByJamm                = sum(sumTable.MissByJamm,'omitnan');

stats.NumOfPulsesAnalyzed = numel(unique(sumTable.pulseNum));
stats.DetectionProb       = stats.TruePositive ./ stats.TotalTargetsDetections;
stats.FalseAlramProb      = stats.nFalseAlarms ./ stats.TotalNonTargets;
stats.MissesDetectionProb = stats.nMissesByFB ./ stats.TotalTargetsDetections;

% % meanyTarget = varfun(@nanmean, sumTable, 'InputVariables', 'corrCoeffFR','GroupingVariables','DetectedAndTrue')
% % boxplot()
figure
for k = 1:numel(ffn)
    subplot(nRows, nCols, k)
    scatter(sumTable.isTarget, sumTable.(ffn{k}), [], sumTable.isClassified, 'filled');
    hold on
    xlim([-0.2, 1.2])
%     plot(sumTable.isTarget(logical(sumTable.FalseAlarm)), sumTable.(ffn{k})(logical(sumTable.FalseAlarm)), 'or' )
    plot(sumTable.isTarget(logical(sumTable.MissByJamm)), sumTable.(ffn{k})(logical(sumTable.MissByJamm)), 'dr' )
    title(ffn{k},'Interpreter','none')
    xlabel('isTarget')
end
 sgtitle('Calssifier results, Color by Clasiification')

 figure
 