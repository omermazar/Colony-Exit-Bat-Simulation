function [Table_FullExperiment ] = CaveExit_Experiment()

% This script tests the effect of the Pd on the Bats' Performance
% Nominal Variables (nBat =5, nPrey =10)
% Input- NumOfTrials : Number of trials in each set of parameters
% Output- save files of each experiment and tsummary table

%% Nominal Parameters load
runPath = 'D:\Omer\BatSimulation\All Directories_backup\Bat Simulation';
cd(runPath)

NominalParams_FilePathName = ...
    'D:\Dropbox\University\מחקר\experiments\caveExit\Code\DefaultParamsTable_CaveExit_Rhinopoma.xlsx';
% load(NominalParams_FilePathName);
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

Params = ReadDefaultParamsTable(NominalParams_FilePathName);

% DefaultParameters is a MAT-file with the parameters neede to run the simulation
% SimParams        = Params.SimParams;
% BatFlightParams  = Params.BatFlight;
% BatSonarParams   = Params.BatSonar;
% TerrainParams    = Params.Terrain;
% PreyFlightParams = Params.PreyFlight;

% Params = ReadDefaultParamsTable(NominalParams_FilePathName);
% The terrain
AllParams.TerrainParams    = Params.Terrain;
AllParams.SimParams        = Params.SimParams;
AllParams.BatFlightParams  = Params.BatFlight;
AllParams.BatSonarParams   = Params.BatSonar;
AllParams.PreyFlightParams = Params.PreyFlight;


%% special params 
%%% General Set-up for the experiment
AllParams.BatSonarParams.TerminalFreqVariance = 0;
AllParams.BatSonarParams.UniqueFreqsFlag      = 'Random'; %'Random'; %'MaxDiff'; %'Random'; % 'MaxDiff'

AllParams.SimParams.TotalBatsNumber       = 1;
AllParams.SimParams.TotalPreysNumber      = 0;
AllParams.PreyFlightParams.React2BatsFlag = 0;
AllParams.BatSonarParams.ReceiverType     = 'FilterBank'; % 'FilterBank'
AllParams.SimParams.TestMode              = 'caveExit';
AllParams.SimParams.DetectObs             = true;
AllParams.SimParams.DetectPrey            = false;
AllParams.SimParams.DetectConsps          = true; %false; %true;
AllParams.SimParams.MaskingByClutter      = false;
AllParams.SimParams.MaskingByConsps       = true;
AllParams.SimParams.MaskingByConspsObsEchoes  = false; %false; % true
AllParams.SimParams.AcousticsCalcFlag     = true;
AllParams.SimParams.SimulationTime        = 15;
AllParams.TerrainParams.Type              = 'Room-L obs'; %'Simple Cave'; %'Room-L obs';
AllParams.BatSonarParams.ClutterPointReflecetingArea = 0.1; %m2


 %% Bulid the Teerrain with the Cave
[Terrain, TerrainParams.Xmax, TerrainParams.Ymax, TerrainParams.Zmax,...
    TerrainParams.Xmin, TerrainParams.Ymin, TerrainParams.Zmin] = ...
    BuildEnvironment(AllParams.TerrainParams, AllParams.SimParams.xyResolution, AllParams.TerrainParams.Type);


%% Parameters modified during the experiment

%%% Var1 - Number Of Bats
% TotalBatsNumber =  [1 2 5 10 15 20 25 40 100]
TotalBatsNumber =  3; %[1 2 5 10 20 40 100] ;% [1 2 5 10 20 40 100]; % [1 2 5 10  20 ]; %; % 20; % 100; %[1 2 5 10 15 20 25 40]; %[ 18 20 22]; %[20 20]; %[1 2 5 10 15 20 25 40];
NumOfSetups1    = numel(TotalBatsNumber); % ([Preg, Post])NumOfSetups1 = 2; % ([Preg, Post])
% % num of repetitions in each setup
nTrials         =  1; %[240 120 48 24 12 6 6];
% nTrials         = 6; % [240 120 48 24 12 6 6];; %[240 120 24 12 6]; % [180 96 36 15 12 8 8 6];% [ 4 4 4]; %[4 4 ]; %   [160 80 32 16 8 4] ; % [40 20 8 4]; %[4 2 1 1]

%%% Var2 - RecieverType and Y/N Masking (if no Acoustics there is no masking
% RecieverType = 'FilterBank' ; % {'CorrelationDetector', 'CorrelationDetector', 'FilterBank'};  {'CorrelationDetector'};
AcousticsMaskFlag = [true] ;  % false; %  ; % [true, false]  ; %[true, false]; % [ 0 1 1];
NumOfSetups2      = length(AcousticsMaskFlag);
AllParams.SimParams.MaskingByConspsConspsEchoes = true;
AllParams.BatSonarParams.TwoEarsReception       = 1;
AllParams.BatSonarParams.ClutterSideEstimate    = 1;

% % Var3
ObsMemorySize = [10, 2]; %[0:4];
NumOfSetups3 = length(ObsMemorySize);


%%% parameter for Saving DATA
SavePath       = 'D:\Dropbox\University\מחקר\experiments\caveExit\Results\';
SaveExpDir     = 'batwoman\test\'; %'Sep2022\Run\17Oct\';
saveBatDATAflg = 0;

% The experiment
ExperimentName = 'Rhino_FalseAlarms_';
% Nexp=0;
Table_FullExperiment = table(); % output table
TableSaveName = ['Table_',ExperimentName];



%%     

% NumOfSetups = length(NumOfBats);
% nFlightNum = Nexp; % Serial Number of the Flight
TotalExperiments = sum(nTrials)*NumOfSetups2; %sum(nTrials) * NumOfSetups1; % * NumOfSetups2;
NumOfCompletedExp = 0;
ntot = 0;

disp([10, 'number of trials: ', num2str(NumOfCompletedExp), ' of: ' , num2str(TotalExperiments)]);
disp([ 'NumberOfBats: ', num2str(AllParams.SimParams.TotalBatsNumber), ...
    ' NumOfPreys: ', num2str(AllParams.SimParams.TotalPreysNumber),...
    ' ReceiverType: ',  num2str(AllParams.BatSonarParams.ReceiverType) ]);
display(datetime(now,'ConvertFrom','datenum'));

for nVar1 = 1: NumOfSetups1 % Change the explanatory varaiables
    AllParams.SimParams.TotalBatsNumber = TotalBatsNumber(nVar1);
    
    for nVar2 = 1:NumOfSetups2 % change the 2nd variablesum(nTrials)
        AllParams.SimParams.MaskingByConsps          = AcousticsMaskFlag(nVar2);
        AllParams.SimParams.MaskingByConspsObsEchoes = AcousticsMaskFlag(nVar2);
        AllParams.SimParams.DetectConsps             = AcousticsMaskFlag(nVar2);

        for nVar3 = 1:NumOfSetups3
            AllParams.BatSonarParams.ObsMemorySize = ObsMemorySize(nVar3);           
            ntot = ntot +1 ;
            
            %%% Set-UP Paramters
              SetupParameters = struct( ...
                'ExperimentName' ,           [ExperimentName, num2str(AllParams.SimParams.TotalBatsNumber), 'Bats_', ...
                                                num2str(AllParams.SimParams.MaskingByConsps), 'Masking_', ...
                                                num2str(AllParams.BatSonarParams.ObsMemorySize), 'Memory_', num2str(ntot) ], ...
                'NumberOfBats',              AllParams.SimParams.TotalBatsNumber, ...
                'NumOfPreys',                AllParams.SimParams.TotalPreysNumber,...
                'TestMode',                  AllParams.SimParams.TestMode , ...
                'ReceiverType',              AllParams.BatSonarParams.ReceiverType, ... % AllParams.BatSonarParams.ReceiverType,...'JAR_Corr'
                'RoomType',                  AllParams.TerrainParams.Type , ...
                'DetectConsps',              AllParams.SimParams.DetectConsps, ...
                'DetectObs',                 AllParams.SimParams.DetectObs, ...
                'DetectPrey',                AllParams.SimParams.DetectPrey, ...
                'MaskingByClutter',          AllParams.SimParams.MaskingByClutter, ...
                'MaskingByConsps',           AllParams.SimParams.MaskingByConsps, ...
                'MaskingByConspsObsEchoes',  AllParams.SimParams.MaskingByConspsObsEchoes, ...
                'MaskingByConspsEchoes',     AllParams.SimParams.MaskingByConspsEchoes, ...
                'MaskingByConspsPreyEchoes', AllParams.SimParams.MaskingByConspsPreyEchoes, ...
                'React2BatsFlag',            AllParams.PreyFlightParams.React2BatsFlag, ...
                'ClutterGainDB',             AllParams.BatSonarParams.ClutterGainDB, ...
                'PulsePower',                AllParams.BatSonarParams.PulsePower, ...
                'ObsMemorySize',             AllParams.BatSonarParams.ObsMemorySize, ...
                'VelocityFree',              AllParams.BatFlightParams.NominalVelocity, ...
                'VelocityManuever',          AllParams.BatFlightParams.MaxVelocity ...
                );
            if  AllParams.BatSonarParams.JARMode
                SetupParameters.ReceiverType = 'JAR_Corr';
            end% if JARMode

            SetupParamFields = fields(SetupParameters);
            %%% Init Results Struct
            SummaryVectorsCell = cell(nTrials(nVar1),1);
            expRawData         = cell(nTrials(nVar1),1);
            BatDATAcell        = cell(nTrials(nVar1),1);
            %%% RUN the Experiment with the required setup
            %PARFOR
            %         nTrials = nTrials(nVar1);


            disp(SetupParameters.ExperimentName);
            display(datetime(now,'ConvertFrom','datenum'));
            % parfor nTrialCounter = 1:nTrials(nVar1)
            for nTrialCounter = 1:nTrials(nVar1)  % Var2NumOfTrials(nVar2) % repeats in each trial

                DataToAnalyze = BatFlightForGui(AllParams, Terrain);

                SummaryVectorsCell{nTrialCounter} = DataToAnalyze.FlightInterferceSummary.SummaryDataVectors;
                if saveBatDATAflg
                    BatDATAcell{nTrialCounter} = DataToAnalyze;
                end
                %             expRawData{nTrialCounter}         = [DataToAnalyze.BAT.InterReportStrctOnLine];
                display(datetime(now,'ConvertFrom','datenum'))
                display(['XXXXX  ', SetupParameters.ExperimentName, ' trial: ', num2str(nTrialCounter),]) 
            end % parfor nTrialCounter = 1:NumOfTrials

            VectorFieldsNames = fields(SummaryVectorsCell{1});

            %%%% extracting the output from the cell
            % Init The Output Struct
            TrialBatsNumber = AllParams.SimParams.TotalBatsNumber;
            NumOfTestedBats = nTrials(nVar1)*TrialBatsNumber; % Var2NumOfTrials(nVar2) * TrialBatsNumber;

            % Insert SetUp Parameters
            for fSetUPIdx = 1:numel(SetupParamFields)
                FName= SetupParamFields{fSetUPIdx};
                % differnt assignment for numeric and strings
                if isnumeric(SetupParameters.(FName)) || isstring(SetupParameters.(FName)) || islogical(SetupParameters.(FName))
                    ExperimetStepSum.(FName)(1:NumOfTestedBats,1) = deal(SetupParameters.(FName));
                elseif  ischar(SetupParameters.(FName))
                    [ExperimetStepSum.(FName){1:NumOfTestedBats,1}] = deal(SetupParameters.(FName));
                end %if isnumeric(SetupParameters.(FName))
            end

            % Insert Result vectors
            for kk = 1:nTrials(nVar1) % Var2NumOfTrials(nVar2)
                for Fidx = 1:numel( VectorFieldsNames)
                    FN = VectorFieldsNames{Fidx};
                    ExperimetStepSum.(FN)( (kk-1)*TrialBatsNumber+1 : kk*TrialBatsNumber, 1)...
                        = SummaryVectorsCell{kk}.(FN)';
                end % for Fidx
            end % for kk

            NumOfCompletedExp = NumOfCompletedExp + nTrials(nVar1); % Var2NumOfTrials(nVar2);
            disp([10, 'number of trials: ', num2str(NumOfCompletedExp), ' of: ' , num2str(TotalExperiments)]);
            disp([ 'NumberOfBats: ', num2str(AllParams.SimParams.TotalBatsNumber), ...
                ' NumOfPreys: ', num2str(AllParams.SimParams.TotalPreysNumber),...
                ' ReceiverType: ',  num2str(AllParams.BatSonarParams.ReceiverType) ]);
            display(datetime(now,'ConvertFrom','datenum'));

            % Concatanated Table of All Trials

            Table_FullExperiment  = [Table_FullExperiment ; struct2table(ExperimetStepSum) ];

            % Save Results
            try
                %             currExpRawData = expRawData;
                save([SavePath, SaveExpDir,SetupParameters.ExperimentName], 'ExperimetStepSum');
                save([SavePath, SaveExpDir,SetupParameters.ExperimentName,'_AllParams'], 'AllParams');
                if saveBatDATAflg
                    save([SavePath, SaveExpDir,SetupParameters.ExperimentName,'_BatDATA'], 'BatDATAcell');
                end %if
                %             save([SavePath, SaveExpDir,SetupParameters.ExperimentName,'_rawData'], 'currExpRawData');
            catch
                hops= 'save ???'
            end
            clear ExperimetStepSum;

        end % for nVar3
    end %for nVar2
        
%     end %for nVar2 = 1:NumOfSetups2 % change the 2nd variable
end % nVar1 = 1:NumOfSetups

% Save Concatenated Tabel
writetable(Table_FullExperiment, [SavePath, SaveExpDir, TableSaveName] );


