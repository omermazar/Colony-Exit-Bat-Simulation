
function [] = myClassifierResultsPlot(FindsMaskingStruct)

% fn = fieldnames(classifierResults{1});
% for n = 1:numel(fn)
%     for k =1:numel(classifierResults)
%         results.(fn{n})(k) = classifierResults{k}.(fn{n});
%     end % for k
% end % for n

classifierResults = FindsMaskingStruct.ClassifierFullResults;

figure; hold on
legend on
% ffn = {'pks_trghs_rms', 'corrCoeffFR', 'corrCoeffMatrix', 'mseFR', 'mseMatrix', 'ssimMarix'};
ffn = { 'corrCoeffFR', 'corrCoeffMatrix'};

for n = 1:numel(ffn)
    if ~strcmp(ffn{n}, 'ntime') && ~strcmp(ffn{n}, 'isClassified')
        plot(classifierResults.ntime, classifierResults.(ffn{n}), 'o-', 'LineWidth', 1.5, DisplayName=ffn{n} ) 
%         pause
    end    
end 
yl = ylim;
plot(FindsMaskingStruct.ClassifierFalseAlarmsTimes, yl(1) * ones(size(FindsMaskingStruct.ClassifierFalseAlarmsTimes)), 'dr', 'MarkerSize', 6, 'LineWidth', 2, 'DisplayName', 'FALSE')
plot(FindsMaskingStruct.ClassifierMaskedTimes,      yl(1) * ones(size(FindsMaskingStruct.ClassifierMaskedTimes)),      '*r', 'MarkerSize', 6, 'LineWidth', 2,'DisplayName', 'Jammed')
plot(FindsMaskingStruct.ClassifierDetecdtedAndTRUE, yl(1) * ones(size(FindsMaskingStruct.ClassifierDetecdtedAndTRUE)), '*g', 'MarkerSize', 6, 'LineWidth', 2, 'DisplayName', 'TRUE')
grid on
title(['Classifier results, Pulse:', num2str(FindsMaskingStruct.PulseNum)])