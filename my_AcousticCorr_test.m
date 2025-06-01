
function [acous_stats] = my_AcousticCorr_test(BatDATA, fields_to_analyze)

% Nov2021

if nargin < 2
    fields_to_analyze = {'DetectedPreys', 'MaskedPreys', 'DetecetedPrey2InterferenceRatioDB', 'SelfCorrealationMaxdB', 'InterCorrelationMaxdB'};
end % if nargin

nfields = numel(fields_to_analyze);

BatNum = 1:BatDATA.AllParams.SimParams.TotalBatsNumber;


acous_stats.BatNum = BatNum';
for ifield = 1:nfields
    curr_field = fields_to_analyze{ifield};
    acous_stats.(strcat(curr_field,'_orig')) = zeros(numel(BatNum),1);
    acous_stats.(strcat(curr_field,'_Acous')) = zeros(numel(BatNum),1);
    
    if strcmp(curr_field, 'DetectedPreys') || strcmp(curr_field, 'MaskedPreys')
        for kBat = BatNum
            acous_stats.(strcat(curr_field,'_orig'))(kBat) = numel( [BatDATA.BAT(kBat).FindsMaskingStruct.(curr_field)] );
            if ~strcmp(curr_field, 'MaskedPreys')
                acous_stats.(strcat(curr_field,'_Acous'))(kBat) = numel( [BatDATA.BAT(kBat).FindsMaskingStruct_Acoustics.(curr_field)] );   
            else
                acous_stats.(strcat(curr_field,'_Acous'))(kBat) = numel( [BatDATA.BAT(kBat).FindsMaskingStruct_Acoustics.(curr_field)] ); 
            end
        end % for kBat
        
    else % if strcmp(curr_field, 'DetectedPreys') || strcmp(curr_field, 'MaskedPreys')
        for kBat = BatNum
            acous_stats.(strcat(curr_field,'_orig'))(kBat) = median( [BatDATA.BAT(kBat).FindsMaskingStruct.(curr_field)] );
            acous_stats.(strcat(curr_field,'_Acous'))(kBat) = median( [BatDATA.BAT(kBat).FindsMaskingStruct_Acoustics.(curr_field)] );     
        end % if strcmp(curr_field, 'DetectedPreys') || strcmp(curr_field, 'Ma
        
    end % for kBat
end % for ifield = 1:nfields

acous_stats.masking_ratio_orig = acous_stats.MaskedPreys_orig ./ acous_stats.DetectedPreys_orig;
acous_stats.masking_ratio_Acous = acous_stats.MaskedPreys_Acous ./ acous_stats.DetectedPreys_Acous;

acous_stats = struct2table(acous_stats);
acous_stats = movevars(acous_stats, 'masking_ratio_orig', 'After', 'MaskedPreys_Acous');
acous_stats = movevars(acous_stats, 'masking_ratio_Acous', 'After', 'masking_ratio_orig');