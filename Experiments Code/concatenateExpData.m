function [tableExperiment, summaryData] = concatenateExpData(pathName, groupBy)

% ConcatenateExpData returns one struct of all the data collected and saved
% for each trials
% Input - PathName- the Directory of the files from the experiment
%           ('D:\Dropbox\University\מחקר\experiments\caveExit\Results\Sep2022\Process')
%         groupBy- A cell with the names of the variables that the
%               suumaryDATA is grouped by, for example: 
%               groupBy = {'NumberOfBats', 'MaskingByConsps', 'ClutterGainDB' , 'ObsMemorySize'}
% Output- ExperimentDATAStruct - a collum- struct with the fields:
            % Inputs - FilghtNum, ProbToDetect,  NumberOfBats, NumOfPreys,  
            % OutPuts - NumOfCatches, NumberOfDetections,NumberOfPulses, NumberOfSuccesulSearchStages
            % DetectionWhileSearching, DetectionWhileApproaching,NumberOfPulsesTryingToHunt

if nargin <2
    groupBy = {'NumberOfBats', 'MaskingByConsps'};
end
            
%%% Grabing the file names in the direcrory
Files = dir(fullfile(pathName,'*.mat'));
fileNames= {Files(:).name};

% New fiels: exist only in part of files
newFields = {'dist2Consps_Mean', 'dist2Consps_Med', 'dist2Consps_Min', 'dist2Consps_STD', 'dist2Consps_STDMed', 'dist2Consps_STDMin'};

% Init
tableExperiment = table();
nFiles = 0;
for k = 1:numel(Files)
    % remove AllPArams
    if ~contains(fileNames{k},'AllParams') && ~contains(fileNames{k},'rawData') && ~contains(fileNames{k},'BatDATA')
        % load file
        currStr = load(strcat(pathName,'\',Files(k).name));
        % add to table
        try
            tt = struct2table(currStr.ExperimetStepSum);
        catch
            warning(strcat('File is not included: ', Files(k).name, ' !!!!'))
        end
        
        %%% check for new field and add them if they are missing
        if ~ismember(newFields{1}, tt.Properties.VariableNames)
            tt{:,newFields} = nan(size(tt,1),numel(newFields));
        end % if ~ismember()
        
        %%% Add Missing fields in partial detection
        missingFields1 = setdiff( tableExperiment.Properties.VariableNames, tt.Properties.VariableNames);
        if ~isempty(missingFields1)
            tt{:,missingFields1} = nan(size(tt,1),numel(missingFields1));
        end
        
        missingFields2 = setdiff( tt.Properties.VariableNames, tableExperiment.Properties.VariableNames);
        if ~isempty(missingFields2)
            tableExperiment{:,missingFields2} = nan(size(tableExperiment,1),numel(missingFields2));
        end

        %%% Concatanate the tables
        try
            tableExperiment = [tableExperiment; tt];
        catch
            warning(['File ',  Files(k).name, ' not Includeded !!!!'])
        end %try

        if mod(nFiles, 10) == 0
            display(datetime(now,'ConvertFrom','datenum','Format','HH:mm:ss'))
            display(['file number ', num2str(nFiles), ' out of ', num2str(numel(Files))])
        end
        
        nFiles = nFiles + 1;
    end % if ~contains()
end % for k

%     g = findgroups(tableExperiment.NumberOfBats);
% [g, IDbats, IDmask, IDdetCons]  = findgroups(tableExperiment.NumberOfBats, tableExperiment.MaskingByConsps, tableExperiment.DetectConsps);
[g, TID]= findgroups(tableExperiment(:, groupBy));

summaryData.totalFiles       = nFiles;
summaryData.fileNames        = sum(~contains (fileNames,'AllParams'));
summaryData.totalExperiments = size(tableExperiment,1);
summaryData.numOfOBats       = unique(tableExperiment.NumberOfBats)';
summaryData.MaskingByConsps  = unique(tableExperiment.MaskingByConsps)';
summaryData.DetectConsps     = unique(tableExperiment.DetectConsps)';

% summaryData.groupsTable     = table(IDbats, IDmask, IDdetCons, splitapply(@numel, tableExperiment.NumberOfBats, g));
% summaryData.groupsTable.Properties.VariableNames{4} = 'total';
summaryData.groupsTable = TID;
summaryData.groupsTable.numel = splitapply(@numel, tableExperiment.NumberOfBats, g);
summaryData.groupsTable.meanExitSuccess  = splitapply(@nanmean, tableExperiment.ExitSuccess, g) ;
summaryData.groupsTable.stdExitSuccess   = splitapply(@nanstd , tableExperiment.ExitSuccess, g) ;