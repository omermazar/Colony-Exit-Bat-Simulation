function [MeanSTDStruct] = CalculateMeanStdForStruct(StructIn, DesiredFieldNames)
% returns the mean and the std of DesiredFieldNames in StructIn

% MeanSTDStruct(length(DesiredFieldNames)) = struct();
AllFieldNames = fieldnames(StructIn);
for k = 1:length(DesiredFieldNames)
    MeanTxt= ['Mean',(DesiredFieldNames{k})];
    StdTxt= ['STD',(DesiredFieldNames{k})];
    if isfield (StructIn, DesiredFieldNames{k})
        MeanSTDStruct.(MeanTxt) = mean(StructIn.(DesiredFieldNames{k}));
        MeanSTDStruct.(StdTxt) = std(StructIn.(DesiredFieldNames{k}));
    else % if isfield (StructIn, DesiredFieldNames{k})
        MeanSTDStruct.(MeanTxt) = [];
        MeanSTDStruct.(StdTxt) = [];
    end % if isfield (StructIn, DesiredFieldNames{k})
end