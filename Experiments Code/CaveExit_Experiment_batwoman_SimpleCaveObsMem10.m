function [Table_FullExperiment ] = CaveExit_Experiment_batwoman_SimpleCaveObsMem10()

% This script tests the effect of the Pd on the Bats' Performance
% Nominal Variables (nBat =5, nPrey =10)
% Input- NumOfTrials : Number of trials in each set of parameters
% Output- save files of each experiment and tsummary table


% nohup matlab -nodisplay -nodesktop -nosplash -nojvm -r "cd('/home/omermazar/caveExit_Experiment'); t=CaveExit_Experiment_batwoman;exit" > log.txt &

%% Nominal Parameters load
% runPath = 'D:\Omer\BatSimulation\All Directories_backup\Bat Simulation';
% NominalParams_FilePathName = ...
%     'D:\Dropbox\University\????????\experiments\caveExit\Code\DefaultParamsTable_caveExit.xlsx';


runPath =  '/home/omermazar/Bat Simulation'; %'D:\Omer\BatSimulation\All Directories_backup\Bat Simulation';

%%%%%% Select which params to run - Important %%%%%%%%%%%%%%%%%%
% NominalParams_FilePathName = '/home/omermazar/caveExit_Experiment/DefaultParamsTable_CaveExit.xlsx';
NominalParams_FilePathName = '/home/omermazar/caveExit_Experiment/DefaultParamsTable_CaveExit_FInal_PK.xlsx';
%%%%%%                                  %%%%%%%%%%%%%%%%%%%%

% warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

addpath(genpath('/home/omermazar/caveExit_Experiment'));
addpath(genpath(runPath));
cd(runPath)
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
AllParams.TerrainParams.Type              = 'Simple Cave'; %'Simple Cave'; %'Room-L obs';
AllParams.BatSonarParams.ClutterPointReflecetingArea = 0.1; %m2

%% Bulid the Terrain with the Cave
[Terrain, TerrainParams.Xmax, TerrainParams.Ymax, TerrainParams.Zmax,...
    TerrainParams.Xmin, TerrainParams.Ymin, TerrainParams.Zmin] = ...
    BuildEnvironment(AllParams.TerrainParams, AllParams.SimParams.xyResolution, AllParams.TerrainParams.Type);


%% Parameters modified during the experiment
getname = @(x) inputname(1);

%%% Var1 - Number Of Bats
% TotalBatsNumber =  [1 2 5 10 15 20 25 40 100]
TotalBatsNumber = [1 2 5 10 20 40 100] ;% [1 2 5 10 20 40 100]; % [1 2 5 10  20 ]; %; % 20; % 100; %[1 2 5 10 15 20 25 40]; %[ 18 20 22]; %[20 20]; %[1 2 5 10 15 20 25 40];
NumOfSetups1    = numel(TotalBatsNumber); % ([Preg, Post])NumOfSetups1 = 2; % ([Preg, Post])
% % num of repetitions in each setup
nTrials         =  [240 120 48 24 12 6 6];
% nTrials         = 6; % [240 120 48 24 12 6 6];; %[240 120 24 12 6]; % [180 96 36 15 12 8 8 6];% [ 4 4 4]; %[4 4 ]; %   [160 80 32 16 8 4] ; % [40 20 8 4]; %[4 2 1 1]
var1 = getname(TotalBatsNumber);

%%% Var2 - RecieverType and Y/N Masking (if no Acoustics there is no masking
% RecieverType = 'FilterBank' ; % {'CorrelationDetector', 'CorrelationDetector', 'FilterBank'};  {'CorrelationDetector'};
AcousticsMaskFlag = true;  % false; %  ; % [true, false]  ; %[true, false]; 
NumOfSetups2      = length(AcousticsMaskFlag);
var2 = getname(AcousticsMaskFlag);


% % Var3
ObsMemorySize = 10 ;% [3, 5, 8, 10]; % default 6
% MaxVelocity     = 5; %[1.5, 2, 4, 5]; % defualt 3
NumOfSetups3    = length(ObsMemorySize);
var3            = getname(ObsMemorySize);

%%% parameter for Saving DATA
SavePath   = '/home/omermazar/caveExit_RelutsXX'; %'D:\Dropbox\University\????????\experiments\caveExit\Results\';
SaveExpDir = '/'; % 'Sep2022\Run\17Oct\';
saveBatDATAflg = 0;

% The experiment
ExperimentName = 'CaveExit_PK_SimpleCave_Mem10';
% Nexp=0;
Table_FullExperiment = table(); % output table
TableSaveName = ['Table_',ExperimentName];


%% Start
disp(' XXXXXXXXX Start Experiment XXXXXXXXXXXXX')
disp(' ')
disp(['Experiment: ', ExperimentName])
display(datetime(now,'ConvertFrom','datenum'));
disp(['Var1: ',  var1,  ': ', num2str(eval(var1)) ]); 
disp(['Repeat: ', num2str(nTrials)]); 
disp(['Var2: ',  var2,  ': ', num2str(eval(var2))]); 
disp(['Var3: ',  var3,  ': ', num2str(eval(var3))]); 
disp('XXXXXXXXXXXXXXXXXXXXXXX')
disp(' ')

%%

% NumOfSetups = length(NumOfBats);
% nFlightNum = Nexp; % Serial Number of the Flight
TotalExperiments = sum(nTrials) * NumOfSetups2 * NumOfSetups3; %sum(nTrials) * NumOfSetups1; % * NumOfSetups2;
NumOfCompletedExp = 0;
ntot = 0;

for nVar1 = 1: NumOfSetups1 % Change the explanatory varaiables /
    AllParams.SimParams.TotalBatsNumber = TotalBatsNumber(nVar1);

    for nVar2 = 1:NumOfSetups2 % change the 2nd variablesum(nTrials)
        %         ntot = ntot +1 ;
        AllParams.SimParams.MaskingByConsps = AcousticsMaskFlag(nVar2);
        AllParams.SimParams.MaskingByConspsObsEchoes = AcousticsMaskFlag(nVar2);
        AllParams.SimParams.DetectConsps = AcousticsMaskFlag(nVar2);
        AllParams.SimParams.MaskingByConspsConspsEchoes = AcousticsMaskFlag(nVar2);

        for nVar3 = 1:NumOfSetups3
            AllParams.BatSonarParams.ObsMemorySize = ObsMemorySize(nVar3);
%             AllParams.BatFlightParams.MaxVelocity     = MaxVelocity(nVar3);
            ntot = ntot +1 ;
            
            %%% Set-UP Paramters
            SetupParameters = struct( ...
                'ExperimentName' ,           [ExperimentName, num2str(AllParams.SimParams.TotalBatsNumber), 'Bats_', ...
                                                num2str(AllParams.SimParams.MaskingByConsps), 'Masking_', ...
                                                num2str(AllParams.BatSonarParams.ObsMemorySize), 'ObsMem', num2str(ntot) ], ...
                'NumberOfBats',                AllParams.SimParams.TotalBatsNumber, ...
                'NumOfPreys',                  AllParams.SimParams.TotalPreysNumber,...
                'TestMode',                    AllParams.SimParams.TestMode , ...
                'ReceiverType',                AllParams.BatSonarParams.ReceiverType, ... % AllParams.BatSonarParams.ReceiverType,...'JAR_Corr'
                'RoomType',                    AllParams.TerrainParams.Type , ...
                'DetectConsps',                AllParams.SimParams.DetectConsps, ...
                'DetectObs',                   AllParams.SimParams.DetectObs, ...
                'DetectPrey',                  AllParams.SimParams.DetectPrey, ...
                'MaskingByClutter',            AllParams.SimParams.MaskingByClutter, ...
                'MaskingByConsps',             AllParams.SimParams.MaskingByConsps, ...
                'MaskingByConspsObsEchoes',    AllParams.SimParams.MaskingByConspsObsEchoes, ...
                'MaskingByConspsEchoes',       AllParams.SimParams.MaskingByConspsEchoes, ...
                'MaskingByConspsPreyEchoes',   AllParams.SimParams.MaskingByConspsPreyEchoes, ...
                'MaskingByConspsConspsEchoes', AllParams.SimParams.MaskingByConspsConspsEchoes, ...
                'React2BatsFlag',              AllParams.PreyFlightParams.React2BatsFlag, ...
                'ClutterGainDB',               AllParams.BatSonarParams.ClutterGainDB, ...
                'PulsePower',                  AllParams.BatSonarParams.PulsePower, ...
                'ObsMemorySize',               AllParams.BatSonarParams.ObsMemorySize, ...
                'VelocityFree',                AllParams.BatFlightParams.NominalVelocity, ...
                'VelocityManuever',            AllParams.BatFlightParams.MaxVelocity, ...
                'ClutterSideEstimate',         AllParams.BatSonarParams.ClutterSideEstimate ...
                );
            if  AllParams.BatSonarParams.JARMode
                SetupParameters.ReceiverType = 'JAR_Corr';
            end% if JARMode

            SetupParamFields = fields(SetupParameters);
            %%% Init Results Struct
            SummaryVectorsCell = cell(nTrials(nVar1),1);
            %         expRawData         = cell(nTrials(nVar1),1);
            BatDATAcell        = cell(nTrials(nVar1),1);

            %%% RUN the Experiment with the required setup
            %PARFOR
            %         nTrials = nTrials(nVar1);

            disp([10, 'number of trials: ', num2str(NumOfCompletedExp), ' of: ' , num2str(TotalExperiments)]);
            disp(SetupParameters.ExperimentName);
%             disp([ 'NumberOfBats: ', num2str(AllParams.SimParams.TotalBatsNumber), ...
%                 ' NumOfPreys: ', num2str(AllParams.SimParams.TotalPreysNumber),...
%                 ' ReceiverType: ',  num2str(AllParams.BatSonarParams.ReceiverType) ]);
            display(datetime(now,'ConvertFrom','datenum'));


            parfor (nTrialCounter = 1:nTrials(nVar1), 12)  % Var2NumOfTrials(nVar2) % repeats in each trial

                DataToAnalyze = BatFlightForGui(AllParams, Terrain);

                SummaryVectorsCell{nTrialCounter} = DataToAnalyze.FlightInterferceSummary.SummaryDataVectors;
                if saveBatDATAflg
                    BatDATAcell{nTrialCounter} = DataToAnalyze;
                end
                display(datetime(now,'ConvertFrom','datenum'))
                display(['XXXXX  ', SetupParameters.ExperimentName, ' trial: ', num2str(nTrialCounter),]) 
                %             expRawData{nTrialCounter}         = [DataToAnalyze.BAT.InterReportStrctOnLine];

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
                    try
                        ExperimetStepSum.(FN)( (kk-1)*TrialBatsNumber+1 : kk*TrialBatsNumber, 1)...
                            = SummaryVectorsCell{kk}.(FN)';
                    catch
                        warning(['Error writning results:  ExperimetStepSum.(FN)', FN])
                    end % try
                end % for Fidx
            end % for kk

            NumOfCompletedExp = NumOfCompletedExp + nTrials(nVar1); % Var2NumOfTrials(nVar2);
            disp([10, 'number of trials: ', num2str(NumOfCompletedExp), ' of: ' , num2str(TotalExperiments)]);
            disp([ 'NumberOfBats: ', num2str(AllParams.SimParams.TotalBatsNumber), ...
                ' MaskingByConsps: ', num2str(AllParams.SimParams.MaskingByConsps),...
                ' DetectConsps: ',  num2str(AllParams.SimParams.DetectConsps) ]);
            display(['saving: ', ExperimentName])
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


