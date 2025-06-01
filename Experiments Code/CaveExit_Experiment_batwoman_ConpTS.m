function [Table_FullExperiment ] = CaveExit_Experiment_batwoman_ConpTS(exeComputer)

% This script tests the effect of the Pd on the Bats' Performance
% Nominal Variables (nBat =5, nPrey =10)
% Input- NumOfTrials : Number of trials in each set of parameters
%       exeComputer - 'local' | 'batwoman'

% Output- save files of each experiment and tsummary table


% nohup matlab -nodisplay -nodesktop -nosplash -nojvm -r "cd('/home/omermazar/caveExit_Experiment'); t=CaveExit_Experiment_batwoman;exit" > log.txt &

% input
if nargin == 0
     exeComputer = 'batwoman'; % 'local'; % 'batwoman';
end

%% Nominal Parameters load
% runPath = 'D:\Omer\BatSimulation\All Directories_backup\Bat Simulation';
% NominalParams_FilePathName = ...
%     'D:\Dropbox\University\????????\experiments\caveExit\Code\DefaultParamsTable_caveExit.xlsx';
switch  exeComputer
    case 'batwoman'
        %% for batwoman run
        disp(' XXXXXXXXX Hello Experiment world XXXXXXXXXXXXX')
        %%%5 most updated version
        % runPath =  '/home/omermazar/Bat Simulation'; %'D:\Omer\BatSimulation\All Directories_backup\Bat Simulation';
        runPath = '/home/omermazar/Bat Simulation_CaveExit_OJT'; %'/home/omermazar/Bat Simulation_Ver08Feb23';
        % NominalParams_FilePathName = '/home/omermazar/caveExit_Experiment/DefaultParamsTable_CaveExit_FInal_PK.xlsx';
        % For Confusion
        % NominalParams_FilePathName = '/home/omermazar/caveExit_Experiment/DefaultParamsTable_CaveExit_FInal_withConfusion_PK.xlsx';
        NominalParams_FilePathName = '/home/omermazar/caveExit_Experiment/DefaultParamsTable_CaveExit_FInal_PK.xlsx';
        addpath(genpath('/home/omermazar/caveExit_Experiment'));
        addpath(genpath(runPath));
        %%%%%% parameter for Saving DATA
        SavePath   = '/home/omermazar/caveExit_Results_OJT';% '/home/omermazar/caveExit_Results'; 
        SaveExpDir = '/'; % 'Sep2022\Run\17Oct\';

        % limit the number of parallel CPUs
        numWorkers = 2; % 24 for batWoamn


    case 'local'
        %% For local run
        % runPath =  'C:\Users\YossiYNB3\Documents\Omer\Bat Simulation';
        % NominalParams_FilePathName = 'C:\Users\YossiYNB3\Documents\University\Research\experiments\caveExit\Code\DefaultParamsTable_CaveExit_FInal_withConfusion_PK.xlsx';
        
        %%%%% XXXXXX %%%%%% Important
        runPath =  'C:\Users\YossiYNB3\Documents\Omer\Bat Simulation_CaveExit_OJT';
        NominalParams_FilePathName = 'C:\Users\YossiYNB3\Documents\University\Research\experiments\caveExit\Code\DefaultParamsTable_CaveExit_FInal_PK.xlsx';
        
        % parameter for Saving DATA
        SavePath   = 'C:\Users\YossiYNB3\Documents\University\Research\experiments\caveExit\Results\forFinal\PK_obsConfusion\Pre-run'; %'D:\Dropbox\University\????????\experiments\caveExit\Results\';
        SaveExpDir = '/'; % 'Sep2022\Run\17Oct\';

        % limit the number of parallel CPUs
        numWorkers = 10; % 24 for batWoamn

end % switch

%%
saveBatDATAflg = 0;

% warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

% addpath(genpath('/home/omermazar/caveExit_Experiment'));
% addpath(genpath(runPath));
cd(runPath)
if ~isfile(NominalParams_FilePathName)
  
    disp('Worng Param file. Please select the xls file again')
    [fName, pathName] = uigetfile(strcat(runPath,'\*.xls'));
    NominalParams_FilePathName = fullfile(pathName, fName);
end

if ~isfolder(strcat([SavePath, SaveExpDir]))
    disp('Worng Save Path. Please select the Folder again')
    SavePath = uigetdir(strcat(runPath));
    SaveExpDir = '/';
end
disp(strcat('Running on: ', runPath))

%% Load Params
Params = ReadDefaultParamsTable(NominalParams_FilePathName);

% % limit the number of parallel CPUs
% numWorkers = 4; % 24 for batWoamn

% The terrain
AllParams.TerrainParams    = Params.Terrain;
AllParams.SimParams        = Params.SimParams;
AllParams.BatFlightParams  = Params.BatFlight;
AllParams.BatSonarParams   = Params.BatSonar;
AllParams.PreyFlightParams = Params.PreyFlight;


%% special params
%%% General Set-up for the experiment
AllParams.BatSonarParams.ConfusionObs  = 0;
AllParams.BatSonarParams.ObsMemorySize = 10; % New June 2024 - seconds XXXX

%%%%% NewJune 2024
AllParams.BatSonarParams.ObsClusterFlag = 0;
AllParams.BatSonarParams.minClusterTH = 2;

% AllParams.BatSonarParams.TerminalFreqVariance = 0;
% AllParams.BatSonarParams.UniqueFreqsFlag      = 'Random'; %'Random'; %'MaxDiff'; %'Random'; % 'MaxDiff'
% AllParams.SimParams.TotalBatsNumber       = 1;
% AllParams.SimParams.TotalPreysNumber      = 0;
% AllParams.SimParams.TestMode              = 'caveExit';
% AllParams.SimParams.DetectObs             = true;
% AllParams.SimParams.DetectPrey            = false;
% AllParams.SimParams.DetectConsps          = true; %false; %true;
% AllParams.SimParams.MaskingByClutter      = false;
% AllParams.SimParams.MaskingByConsps       = true;
% AllParams.TerrainParams.Type              = 'Room-L obs'; %'Simple Cave'; %'Room-L obs';
% AllParams.BatSonarParams.ClutterPointReflecetingArea = 0.1; %m2

%% Bulid the Terrain with the Cave
[Terrain, TerrainParams.Xmax, TerrainParams.Ymax, TerrainParams.Zmax,...
    TerrainParams.Xmin, TerrainParams.Ymin, TerrainParams.Zmin] = ...
    BuildEnvironment(AllParams.TerrainParams, AllParams.SimParams.xyResolution, AllParams.TerrainParams.Type);


%% Parameters modified during the experiment
getname = @(x) inputname(1);

%%% Var1 - Number Of Bats 
% TotalBatsNumber =  [1 2 5 10 15 20 25 40 100]
TotalBatsNumber = [100] ; % [1 2 5 10 20 40 100]; % [1 2 5 10 ]; % [1 2 5 10 20 40 100]; % [1 2 5 10]; % [1 2 5 10 20 40 100]; 
NumOfSetups1    = numel(TotalBatsNumber); % ([Preg, Post])NumOfSetups1 = 2; % ([Preg, Post])
% % num of repetitions in each setup
nTrials         = [ 2]; % [240 120 48 24 12 6 6]; % [4 4 4 4]; % [ 12 6 6]; [240 120 48 24 12 6 6];   
% nTrials         = 6; % [240 120 48 24 12 6 6];; %[240 120 24 12 6]; % [180 96 36 15 12 8 8 6];% [ 4 4 4]; %[4 4 ]; %   [160 80 32 16 8 4] ; % [40 20 8 4]; %[4 2 1 1]
var1 = getname(TotalBatsNumber);
val1 = eval(var1);

%%% Var2 - Confusion Flag
% only cuncusion and cluster run Important XXXXX
ConfusionObs   = [true, true, false]; % [true, true, false]; % [true, false];  % false; %  ; % [true, false]  ; %[true, false]; % [ 0 1 1];
NumOfSetups2   = length(ConfusionObs);
var2 = getname(ConfusionObs);
val2 = eval(var2);
% Special Case - if CofusionObs is true- test with and without clustering
ClusterObs    = [true, false, false] ; % [true, false, false];  
ObsMemorySize = [1, 10, 10]; % [1, 10, 10]; % with clustering the ObsMemorySize is 1 sec, without it, it is 10 calls

% % Var3 - only for option
TotalPreysNumber = 0 ;
% AllParams.BatSonarParams.ClutterSideEstimate = 1;
% AllParams.BatSonarParams.TwoEarsReception    = 1 ;
% nominalVel = [3, 5, 8, 10];
% manVel = nominalVel ./ 2;
var3 = getname(TotalPreysNumber);
val3 = eval(var3);
NumOfSetups3 = length(val3);



% The experiment
ExperimentName = 'CaveExit_PK_OJT_ConfWithClusterRun100_'; % 'CaveExit_PK_TwoEars_withFA_';
% Nexp=0;
Table_FullExperiment = table(); % output table
TableSaveName = ['Table_',ExperimentName];


%% Start
disp(' XXXXXXXXX Start Experiment XXXXXXXXXXXXX')
disp(' ')
disp(['Experiment: ', ExperimentName])
disp(['runPath: ', runPath])
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
    v1 = val1(nVar1);

    for nVar2 = 1:NumOfSetups2 % change the 2nd variablesum(nTrials)
        %         ntot = ntot +1 ;
       AllParams.BatSonarParams.ConfusionObs = val2(nVar2);
       v2 = val2(nVar2);
       %%% Important - special for this test
       AllParams.BatSonarParams.ObsClusterFlag = ClusterObs(nVar2);
       AllParams.BatSonarParams.ObsMemorySize  = ObsMemorySize(nVar2);
       %%%%
       
       for nVar3 = 1:NumOfSetups3
            % AllParams.BatSonarParams.TwoEarsReception = TwoEarsReception(nVar3);
%             AllParams.BatFlightParams.NominalVelocity = nominalVel(nVar3);
%             AllParams.BatFlightParams.MaxVelocity     = manVel(nVar3);
             
            v3 = val3(nVar3);
            
            ntot = ntot +1 ;
            
            %%% Set-UP Paramters
            SetupParameters = struct( ... 
                'ExperimentName' ,           [ExperimentName, var1, num2str(v1), ... 
                                                '_', var2, num2str(v2),  ...
                                                '_', 'ClusterObs', char(ClusterObs(nVar2)+'0'),  ... %%%%% Cahnged XXXXXXX
                                                '_', num2str(ntot) ], ...
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
                'ConfusionObs',                AllParams.BatSonarParams.ConfusionObs, ...
                'ObsClusterFlag',              AllParams.BatSonarParams.ObsClusterFlag, ...
                'minClusterTH' ,               AllParams.BatSonarParams.minClusterTH, ...
                'VelocityFree',                AllParams.BatFlightParams.NominalVelocity, ...
                'VelocityManuever',            AllParams.BatFlightParams.MaxVelocity, ...
                'ClutterSideEstimate',         AllParams.BatSonarParams.ClutterSideEstimate, ...
                'TwoEarsReception',            AllParams.BatSonarParams.TwoEarsReception ...
                );

            % % New
            % SetupParameters = struct( ...
            %     'ExperimentName' ,           [ExperimentName, var1, num2str(v1), ...
            %                                     '_', var2, num2str(v2),  ...
            %                                     '_', var3, num2str(v3),  ...
            %                                     '_', num2str(ntot) ], ...
            %     'NumberOfBats',                AllParams.SimParams.TotalBatsNumber, ...
            %     'NumOfPreys',                  AllParams.SimParams.TotalPreysNumber,...
            %     'TestMode',                    AllParams.SimParams.TestMode , ...
            %     'ReceiverType',                AllParams.BatSonarParams.ReceiverType, ... % AllParams.BatSonarParams.ReceiverType,...'JAR_Corr'
            %     'RoomType',                    AllParams.TerrainParams.Type , ...
            %     'DetectConsps',                AllParams.SimParams.DetectConsps, ...
            %     'DetectObs',                   AllParams.SimParams.DetectObs, ...
            %     'DetectPrey',                  AllParams.SimParams.DetectPrey, ...
            %     'MaskingByClutter',            AllParams.SimParams.MaskingByClutter, ...
            %     'MaskingByConsps',             AllParams.SimParams.MaskingByConsps, ...
            %     'MaskingByConspsObsEchoes',    AllParams.SimParams.MaskingByConspsObsEchoes, ...
            %     'MaskingByConspsEchoes',       AllParams.SimParams.MaskingByConspsEchoes, ...
            %     'MaskingByConspsPreyEchoes',   AllParams.SimParams.MaskingByConspsPreyEchoes, ...
            %     'MaskingByConspsConspsEchoes', AllParams.SimParams.MaskingByConspsConspsEchoes, ...
            %     'React2BatsFlag',              AllParams.PreyFlightParams.React2BatsFlag, ...
            %     'ClutterGainDB',               AllParams.BatSonarParams.ClutterGainDB, ...
            %     'PulsePower',                  AllParams.BatSonarParams.PulsePower, ...
            %     'ObsMemorySize',               AllParams.BatSonarParams.ObsMemorySize, ...
            %     'VelocityFree',                AllParams.BatFlightParams.NominalVelocity, ...
            %     'VelocityManuever',            AllParams.BatFlightParams.MaxVelocity, ...
            %     'ClutterSideEstimate',         AllParams.BatSonarParams.ClutterSideEstimate, ...
            %     'TwoEarsReception',            AllParams.BatSonarParams.TwoEarsReception, ...
            %     'ConfusionObs',                AllParams.BatSonarParams.ConfusionObs, ...
            %     'ObsClusterFlag',              AllParams.BatSonarParams.ObsClusterFlag, ...
            %     'minClusterTH' ,               AllParams.BatSonarParams.minClusterTH ...
            %     );

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

            disp(SetupParameters.ExperimentName);
%             disp([ 'NumberOfBats: ', num2str(AllParams.SimParams.TotalBatsNumber), ...
%                 ' NumOfPreys: ', num2str(AllParams.SimParams.TotalPreysNumber),...
%                 ' ReceiverType: ',  num2str(AllParams.BatSonarParams.ReceiverType) ]);
            disp(datetime("now"));
            if isempty(gcp('nocreate'))
                parpool(numWorkers);
            end
            parfor (nTrialCounter = 1:nTrials(nVar1), numWorkers) % (nTrialCounter = 1:nTrials(nVar1), numWorkers) % Var2NumOfTrials(nVar2) % repeats in each trial

                DataToAnalyze = BatFlightForGui(AllParams, Terrain);

                SummaryVectorsCell{nTrialCounter} = DataToAnalyze.FlightInterferceSummary.SummaryDataVectors;
                if saveBatDATAflg
                    BatDATAcell{nTrialCounter} = DataToAnalyze;
                end
                display(datetime("now"))
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
            disp(['Finished -  number of trials: ', num2str(NumOfCompletedExp), ' of: ' , num2str(TotalExperiments)]);
            disp([ 'Bats: ', num2str(AllParams.SimParams.TotalBatsNumber), ...
                ' MaskingByConsps: ', num2str(AllParams.SimParams.MaskingByConsps),...
                ' DetectConsps: ',  num2str(AllParams.SimParams.DetectConsps) ]);
            display(['saving: ', ExperimentName])
            display(datetime(now,'ConvertFrom','datenum'));

            % Concatanated Table of All Trials
            try
                Table_FullExperiment  = [Table_FullExperiment ; struct2table(ExperimetStepSum) ];
            catch
                 warning(strcat('Error in adding Data to Table...', var1, '=', num2str(AllParams.SimParams.TotalBatsNumber), ' DAtA is not saved!!!'))
            end
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
                warning('Data is not saved!!!')
                disp('check save path ?')
            end
            clear ExperimetStepSum;

        end % for nVar3

    end %for nVar2

    %     end %for nVar2 = 1:NumOfSetups2 % change the 2nd variable
end % nVar1 = 1:NumOfSetups

% Save Concatenated Tabel
writetable(Table_FullExperiment, [SavePath, SaveExpDir, TableSaveName] );


