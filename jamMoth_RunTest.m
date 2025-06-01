function [Table_FullExperiment ] = jamMoth_RunTest()

% This script tests the effect of the Pd on the Bats' Performance
% Nominal Variables (nBat =5, nPrey =10)
% Input- NumOfTrials : Number of trials in each set of parameters
% Output- save files of each experiment and tsummary table

%% Nominal Parameters load
NominalParams_FilePathName = ...
    'D:\Omer\BatSimulation\All Directories_backup\Bat Simulation_Jamming_Moth\DATA\DefaultParamsTable.xlsx';
% load(NominalParams_FilePathName);
Params = ReadDefaultParamsTable(NominalParams_FilePathName);
% DefaultParameters is a MAT-file with the parameters neede to run the
% simulation
AllParams.SimParams        = Params.SimParams;
AllParams.BatFlightParams  = Params.BatFlight;
AllParams.BatSonarParams   = Params.BatSonar;
AllParams.TerrainParams    = Params.Terrain;
AllParams.PreyFlightParams = Params.PreyFlight;

% The terrain
[Terrain, TerrainParams.Xmax, TerrainParams.Ymax, TerrainParams.Zmax,...
    TerrainParams.Xmin, TerrainParams.Ymin, TerrainParams.Zmin] = ...
    BuildEnvironment(AllParams.TerrainParams, AllParams.SimParams.xyResolution, AllParams.TerrainParams.Type);



%% special params for clutter

%%% General Set-up for the experiment
AllParams.BatSonarParams.ReceiverType     = 'FilterBank';
AllParams.PreyFlightParams.React2BatsFlag = 1;
AllParams.SimParams.TotalPreysNumber      = 5;
AllParams.SimParams.TotalBatsNumber       = 1;
AllParams.SimParams.SimulationTime        = 10;


%% Parameters modified during the experiment

% Var1
JamTxLvl     =   [70 80 85 90 95 100]; %[0 1 1 1 1 1]; %[0 1 1 1 1 1];
NumOfSetups1 = length(JamTxLvl);

Var2NumOfTrials = 20 ; % 20 ; %[ 30 20 12 9 6 3] ; % [ 20 10 6 4 2 2];

% % Var3
% JARMode = 0; % [0 1];
% NumOfSetups3 = length(JARMode);

%%% parameter for Saving DATA
SavePath   = 'D:\Dropbox\University\מחקר\experiments\Moth_Jamming\Results\';
SaveExpDir = 'jamMoth_Sim\';

% The experiment
ExperimentName = 'jamLvl_ShortActive  ';
% Nexp=0;
Table_FullExperiment = table(); % output table
TableSaveName        = ['Table_',ExperimentName];




% NumOfSetups = length(NumOfBats);
% nFlightNum = Nexp; % Serial Number of the Flight
TotalExperiments = sum(Var2NumOfTrials)* NumOfSetups1; % * NumOfSetups2;
NumOfCompletedExp = 0;
ntot = 0;
display("start time:")
display( datetime(now, "ConvertFrom", "datenum") )
for nVar1 = 1:NumOfSetups1 % Change the explanatory varaiables
    tic
    AllParams.PreyFlightParams.JamTxLvl = JamTxLvl(nVar1);

%     for nVar2 = 1:NumOfSetups2 % change the 2nd variable
%         %  ntot = ntot +1 ;
% 
%         AllParams.SimParams.TotalBatsNumber = NumOfBats(nVar2);
%         
%         for nVar2 = 1:NumOfSetups2
            ntot = ntot +1 ;
            
            %%% Set-UP Paramters
            SetupParameters = struct( ...
                'ExperimentName' ,      [ExperimentName, num2str(ntot) ], ...
                'NumberOfBats',         AllParams.SimParams.TotalBatsNumber, ...
                'NumOfPreys',           AllParams.SimParams.TotalPreysNumber,...
                'PreyTOBatsRatio',      AllParams.SimParams.TotalPreysNumber ./  AllParams.SimParams.TotalBatsNumber,...
                'ReceiverType',         AllParams.BatSonarParams.ReceiverType, ... % AllParams.BatSonarParams.ReceiverType,...'JAR_Corr'
                'PulseDuration_Search', AllParams.BatSonarParams.PulseDuration_Search,...
                'PulsePower',           AllParams.BatSonarParams.PulsePower, ...
                'JamTxLvl',             AllParams.PreyFlightParams.JamTxLvl, ...
                'JamSigDetectionTH',    AllParams.PreyFlightParams.JamSigDetectionTH, ...
                'JamActiveStartFreq',   AllParams.PreyFlightParams.JamActiveStartFreq, ...
                'JamActiveEndFreq',     AllParams.PreyFlightParams.JamActiveEndFreq, ...
                'JamActiveDur',     AllParams.PreyFlightParams.JamActiveDur ...
                );
            if  AllParams.BatSonarParams.JARMode
                SetupParameters.ReceiverType = 'JAR_Corr';
            end% if JARMode
            
            SetupParamFields = fields(SetupParameters);
            %%% Init Results Struct
            SummaryVectorsCell = cell(Var2NumOfTrials,1);
            
            %%% RUN the Experiment with the required setup
            %PARFOR
            %% the jamming signal
            AllParams.PreyFlightParams.jamPreySignal = generateAcousicPreyJamming(AllParams.PreyFlightParams, AllParams.SimParams.FsAcoustic);
            AllParams.PreyFlightParams.JamLoadFlg    = 1;

            parfor nTrialCounter = 1:Var2NumOfTrials  % Var2NumOfTrials(nVar2) % repeats in each trial
                
                DataToAnalyze = BatFlightForGui(AllParams, Terrain);
                
                SummaryVectorsCell{nTrialCounter} = DataToAnalyze.FlightInterferceSummary.SummaryDataVectors;
            end % parfor nTrialCounter = 1:NumOfTrials
            
            VectorFieldsNames = fields(SummaryVectorsCell{1});
            
            %%%% extracting the output from the cell
            % Init The Output Struct
            TrialBatsNumber = AllParams.SimParams.TotalBatsNumber;
            NumOfTestedBats = Var2NumOfTrials*TrialBatsNumber; % Var2NumOfTrials(nVar2) * TrialBatsNumber;
            
            % Insert SetUp Parameters
            for fSetUPIdx = 1:numel(SetupParamFields)
                FName= SetupParamFields{fSetUPIdx};
                % differnt assignment for numeric and strings
                if isnumeric(SetupParameters.(FName))
                    ExperimetStepSum.(FName)(1:NumOfTestedBats,1) = deal(SetupParameters.(FName));
                elseif  ischar(SetupParameters.(FName))
                    [ExperimetStepSum.(FName){1:NumOfTestedBats,1}] = deal(SetupParameters.(FName));
                end %if isnumeric(SetupParameters.(FName))
            end
            
            % Insert Result vectors
            for kk = 1:Var2NumOfTrials % Var2NumOfTrials(nVar2)
                for Fidx = 1:numel( VectorFieldsNames)
                    FN = VectorFieldsNames{Fidx};
                    ExperimetStepSum.(FN)( (kk-1)*TrialBatsNumber+1 : kk*TrialBatsNumber, 1)...
                        = SummaryVectorsCell{kk}.(FN)';
                end % for Fidx
            end % for kk
            
            NumOfCompletedExp = NumOfCompletedExp + Var2NumOfTrials; % Var2NumOfTrials(nVar2);
            disp([10, 'number of trials: ', num2str(NumOfCompletedExp), ' of: ' , num2str(TotalExperiments)]);
            disp([ 'NumberOfBats: ', num2str(AllParams.SimParams.TotalBatsNumber), ...
                ' NumOfPreys: ', num2str(AllParams.SimParams.TotalPreysNumber),...
                ' Jam  Lvl: ',  num2str(AllParams.PreyFlightParams.JamTxLvl) ]);
            display( datetime(now, "ConvertFrom", "datenum") )
            % Concatanated Table of All Trials
            
            Table_FullExperiment  = [Table_FullExperiment ; struct2table(ExperimetStepSum) ];
            
            % Save Results
            try
                save([SavePath, SaveExpDir,SetupParameters.ExperimentName], 'ExperimetStepSum');
                save([SavePath, SaveExpDir,SetupParameters.ExperimentName, '_AllParams'], 'AllParams');
            catch
                hops= 'save ???'
            end
            clear ExperimetStepSum;
            toc
%         end %for nVar2
        
%     end %for nVar2 = 1:NumOfSetups2 % change the 2nd variable
end % nVar1 = 1:NumOfSetups

% Save Concatenated Tabel
writetable(Table_FullExperiment, [SavePath, SaveExpDir, TableSaveName] );


