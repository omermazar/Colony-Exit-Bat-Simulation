function [] = MyBatAnimation(DataToAnalyze, DataToPlot, SaveVideoFlag, Terrain, AnimationTimeRate, AnimmateCurrentOnly, DetectedObj, PathName)


% DataToAnalyze  = BatDATA;
% DataToPlot - Struct:
% DataToPlot.AnimateStartTime = 0;
% DataToPlot.AnimateEndTime   = inf;
% DataToPlot.PlotAllBatsFlag  = 1;
% DataToPlot.BatNumToPlot     = 4;
% DataToPlot.PlotAllPreysFlag = 0;
% DataToPlot.PreyNumToPlot    = 0;
% DataToPlot.JamPreyPosFlag  = 1;

% Terrain = BatDATA.AllParams.TerrainParams.terrainPoints;
% AnimationTimeRate = 75; The bigger the faster
% AnimmateCurrentOnly = false; % true/ false
% DetectedObj = 'Prey'; % 'Prey'; 'Obs'; 'Conspecific'

%% input
if nargin < 6
     DetectedObj = 'Prey'; % 'Prey'; 'Obs'; 'Conspecific'
end

if nargin < 5
    AnimationTimeRate = 70;
end % nargin

% DataToPlot - Struct:
% DataToPlot.AnimateStartTime = 0;
% DataToPlot.AnimateEndTime   = inf;
% DataToPlot.PlotAllBatsFlag  = 1;
% DataToPlot.BatNumToPlot     = 1;
% DataToPlot.PlotAllPreysFlag = 0;
% DataToPlot.PreyNumToPlot    = 0;
% DataToPlot.JamPreyPosFlag  = 1;


% save to Video
% SaveVideoFlag = 0;
if SaveVideoFlag
    if ~exist('PathName', 'var')
        %       PathName = 'D:\Dropbox\University\מחקר\experiments\Moth_Jamming\Animation\';
        %       VideoFileName = ['BatVideo', num2str(round(10000*rand(1))), '.avi'];
        PathName = strcat(cd, '\animationMovies\');
    end
    VideoFileName = strcat('BatVideo_', string(datetime(now,'ConvertFrom','datenum','Format','ddMMMuu_HHmmss')), '.avi');

    VideoObj = VideoWriter(strcat(PathName, VideoFileName));
    VideoObj.FrameRate = 180;
    open(VideoObj);
end % if SaveVideoFlag


% Consts
SimParams = DataToAnalyze.AllParams.SimParams;

NumberOfBats = SimParams.TotalBatsNumber;
NumberOfPreys= SimParams.TotalPreysNumber;

SampleTime =  DataToAnalyze.AllParams.SimParams.SampleTime;
MaxTime =  DataToAnalyze.AllParams.SimParams.SimulationTime;
if strcmp(SimParams.TestMode, 'caveExit')
    mm = [DataToAnalyze.BAT.ExitnTime]*SampleTime; 
    mm(isnan(mm)) = MaxTime;
    MaxTime = max(mm)*1.1;
end % if strcmp

TimeVec = [0:SampleTime:MaxTime]; % Time Vector
xyResolution = DataToAnalyze.AllParams.SimParams.xyResolution;

CatchPreyTimesMat = zeros(NumberOfBats, NumberOfPreys);

FsAcoustic = DataToAnalyze.AllParams.SimParams.FsAcoustic;

startXyIx = max([1, round(DataToPlot.AnimateStartTime / SampleTime)]) ;

endXyIx   = min([MaxTime/SampleTime, round(DataToPlot.AnimateEndTime / SampleTime)]) ;

%%% Setting Colors %%%
[colStruct] = myColorScheme(SimParams.TestMode, NumberOfBats, NumberOfPreys);

ColorOfBatPos   = colStruct.BatPos;
ColorOfPreysPos = colStruct.PreysPos;

% ColorForagingMan = [0 0 0]; %Black
ColorApproach    = colStruct.Approach; % [ 0 1 0];  % Green
ColorBuzzMan     = colStruct.BuzzMan; % [ 0 0 0 ]; % black
ColorObsMan      = colStruct.ObsMan; % [0 0.7 0]; % Green
ColorCaveExit    = colStruct.CaveExit; % [0 1 1]; %
ColorAvoidBatMan = colStruct.AvoidBatMan; % 'y';

ColorFindingPrey = colStruct.preyFinding; % [1 0 1]; % magenta
ColorEstimatePrey= colStruct.preyEstimate; % [1 0 1]; % magenta

%% The figure
fig = figure;
fig.Position = [ 300   150   560   600];
h1 = subplot(2,1,1); hold on;

xMin = min(Terrain(1:2:end, :),[],'all')*xyResolution;
xMax = max(Terrain(1:2:end, :),[],'all')*xyResolution;
yMin = min(Terrain(2:2:end, :),[],'all')*xyResolution;
yMax = max(Terrain(2:2:end, :),[],'all')*xyResolution;
if xMin < -15
    xMin = inf;
    xMax = -inf;
    yMin = inf;
    yMax = -inf;
end

% The maximal position and minimal position of the bats
xBatAll = [DataToAnalyze.BAT.xBatPos];
yBatAll = [DataToAnalyze.BAT.yBatPos];
xMin = min(xMin, min(xBatAll));
xMax = max(xMax, max(xBatAll));
yMin = min(yMin, min(yBatAll));
yMax = max(yMax, max(yBatAll));

xlim([xMin, xMax]);
ylim([yMin, yMax]);
if strcmp(DataToAnalyze.AllParams.TerrainParams.Type,'Room-L obs')
    xlim([0, 10]);
    ylim([0, 10]);
end
xlabel(h1, 'X (m)'); ylabel(h1, 'Y (m)')

h2 = subplot(2,1,2); hold on; 
xlim(h2, [0, SimParams.SimulationTime]);
ylim(h2, [0,140]);
xlabel(h2, 'Time (sec)');
ylabel(h2, 'Signal lvl (dB)');
grid(h2, 'on')

h1.Position = [0.0800    0.36    0.8250    0.6 ]; % [0.1300    0.46    0.7750    0.5 ];
h2.Position = [0.0800    0.09    0.8250    0.17 ]; % [0.1300    0.09    0.7750    0.25 ]

%%%The environment
sz = size(Terrain,1);
for i=1:+2:sz-1
    plot(h1, Terrain(i,:)*xyResolution,Terrain(i+1,:)*xyResolution,'b-')
end



%%% Which Bats to plot %%%
if DataToPlot.PlotAllBatsFlag 
    NumberOfBatsToPlot = NumberOfBats;
    BatNumToPlotVec = 1:NumberOfBatsToPlot;
    title(h1, [num2str(NumberOfBats), ' BATS'])
        % title(h1, [num2str(NumberOfBats), ' BATS, ', DetectedObj])

else %if DataToPlot.PlotAllBatsFlag
    NumberOfBatsToPlot = 1;
    BatNumToPlotVec = [DataToPlot.BatNumToPlot];
    title(h1, ['BAT: ', num2str(BatNumToPlotVec)]);
end % if DataToPlot.PlotAllBatsFlag
title(h2, ['BAT: ', num2str(DataToPlot.BatNumToPlot)]);

%%% Which Preys to plot %%%
if DataToPlot.PlotAllPreysFlag
    NumberOfPreysToPlot = DataToAnalyze.AllParams.SimParams.TotalPreysNumber;
    PreyNumToPlotVec = 1:NumberOfPreysToPlot;

else %if DataToPlot.PlotAllBatsFlag

    if DataToPlot.PreyNumToPlot ~= 0
        NumberOfPreysToPlot = 1;
        PreyNumToPlotVec    = [DataToPlot.PreyNumToPlot];
    else % if PreyNumToPlot ~= 0
        NumberOfPreysToPlot = 0;
    end % if PreyNumToPlot ~= 0
end % if DataToPlot.PlotAllPreysFlag
NumberOfPreysToPlot = min(NumberOfPreysToPlot, DataToAnalyze.AllParams.SimParams.TotalPreysNumber);


%% The figure  - Init
% hAnimationFigure = figure();
% axis([ -1,DataToAnalyze.AllParams.TerrainParams.Xmax+1,...
%         -1, DataToAnalyze.AllParams.TerrainParams.Ymax+1]);   

%%%% Pause button
% btn = uibutton(fig,'push',...
%                'Position',[420, 218, 100, 22],...
%                'ButtonPushedFcn', @(btn,event) plotButtonPushed(btn,ax));
   

%%% init the lines and collecting DATA XY %%%
hBatPos{NumberOfBatsToPlot}        = animatedline(h1);
if NumberOfPreysToPlot > 0
    hPreyPos{NumberOfPreysToPlot}  = animatedline(h1);
    hPreyJamPos{NumberOfPreysToPlot}  = animatedline(h1);
end % if NumberOfPreysToPlot > 0
hCatchPreyPos{NumberOfBatsToPlot}     = animatedline(h1);
hApprocahStage{NumberOfBatsToPlot}    = animatedline(h1);
hBuzzStage{NumberOfBatsToPlot}        = animatedline(h1);
hObsManStage{NumberOfBatsToPlot}      = animatedline(h1);
hJammedPos{NumberOfBatsToPlot}        = animatedline(h1);
hCrushesPos{NumberOfBatsToPlot}       = animatedline(h1);
hClutterFinds{NumberOfBatsToPlot}     = animatedline(h1);
hClutterEstimates{NumberOfBatsToPlot} = animatedline(h1);
hAVoitBatManStage{NumberOfBatsToPlot} = animatedline(h1);

hCaveExitStage{NumberOfBatsToPlot}     = animatedline(h1);
hPreyFinds{NumberOfBatsToPlot}         = animatedline(h1);
hPreyEstimates{NumberOfBatsToPlot}     = animatedline(h1);
hConspsFinds{NumberOfBatsToPlot}       = animatedline(h1);


for nn = 1:NumberOfPreysToPlot
    hPreyPos{nn}    = animatedline(h1,'MarkerFaceColor', ColorOfPreysPos(nn,:), 'Marker', 'o', 'MarkerSize', 4,  'MarkerEdgeColor', 'none');
    hPreyJamPos{nn} = animatedline(h1,'MarkerFaceColor', ColorOfPreysPos(nn,:), 'Marker', '+', 'MarkerSize', 9,  'MarkerEdgeColor', 'r', ...
        'LineWidth', 2);
%     hPreyPos{nn} = animatedline(h1,'color', ColorOfPreysPos(nn,:), 'LineWidth',4);
end

hSelectedBatPos = animatedline(h1, 'MarkerFaceColor','b', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'none');
uistack(hSelectedBatPos,'top');

for nn = 1:NumberOfBatsToPlot
    %% XY AXIS
    % the lines
    hBatPos{nn} = animatedline(h1, 'MarkerFaceColor', ColorOfBatPos(nn,:), 'Marker', 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'none');
    
%     hBatPos{nn} = animatedline(h1, 'color', ColorOfBatPos(nn,:), 'LineWidth',8);
    
    hApprocahStage{nn} = animatedline(h1,'Marker', '.', ...
        'color', ColorApproach ,'linewidth',3,'LineStyle','none'); 
    
    hBuzzStage{nn} = animatedline(h1,'Marker', '.', ...
        'color', ColorBuzzMan ,'linewidth',6,'LineStyle','none'); 
    
    hObsManStage{nn} = animatedline(h1,'Marker', '.', ...
        'color', ColorObsMan ,'linewidth',6,'LineStyle','none');
    
     hAVoitBatManStage{nn} = animatedline(h1,'Marker', '.', ...
        'color', ColorAvoidBatMan ,'linewidth',6,'LineStyle','none');

    hCaveExitStage{nn} = animatedline(h1,'Marker', '.', ...
        'color', ColorCaveExit ,'linewidth',3,'LineStyle','none'); 
    
    hJammedPos{nn} = animatedline(h1,'Marker','d',...
        'MarkerSize',8,'color',colStruct.masking, 'linewidth',2,'LineStyle','none');
    
    hCatchPreyPos{nn} = animatedline(h1,'color', colStruct.catchPrey, 'Marker','+',...
        'MarkerSize',12, 'linewidth',4, 'LineStyle','none'); 
    
    hCrushesPos{nn} = animatedline(h1,'Marker','x',...
        'MarkerSize',12,'color', colStruct.crushes,'linewidth',3,'LineStyle','none');
    
    hClutterFinds{nn}  = animatedline(h1, 'color', colStruct.obsFinding, 'Marker','o',...
        'MarkerSize', 5, 'linewidth',2, 'LineStyle','none');
    
    hClutterEstimates{nn} = animatedline(h1, 'MarkerFaceColor', colStruct.obsEstimate, ...
        'MarkerEdgeColor', colStruct.obsEstimate, 'Marker','s',...
        'MarkerSize', 6, 'linewidth',2, 'LineStyle','none');

    hPreyFinds{nn}  = animatedline(h1, 'color', colStruct.preyFinding, 'Marker','o',...
        'MarkerSize', 5, 'linewidth',1, 'LineStyle','none');
    
    hPreyEstimates{nn} = animatedline(h1, 'MarkerFaceColor', colStruct.preyEstimate, ...
        'MarkerEdgeColor', colStruct.preyEstimate, 'Marker','s',...
        'MarkerSize', 6, 'linewidth',1, 'LineStyle', 'none');

    hConspsFinds{nn} = animatedline(h1, 'MarkerFaceColor', colStruct.conspEstimate, ...
        'MarkerEdgeColor', colStruct.conspEstimate, 'Marker','*',...
        'MarkerSize', 12, 'linewidth',2, 'LineStyle', 'none');
    
    %% Time AXIS
    if SimParams.AcousticsCalcFlag
        hCallsTimePlot{nn}          = animatedline(h2, 'Color', colStruct.ownCalls, 'LineWidth', 2, ...
            'MaximumNumPoints', max(numel(DataToAnalyze.BAT(1).AcousticSig_All),1) ); % calls
        hWantedEchoesTimePlot{nn}   = animatedline(h2, 'Color', colStruct.wantedSignals,  'LineWidth', 2, ...
            'MaximumNumPoints', max(numel(DataToAnalyze.BAT(1).AcousticSig_All),1)); % wanted echoes echoes
        hMaskingSignalsTimePlot{nn} = animatedline(h2, 'Color', colStruct.conspsCalls,  'LineWidth', 2, ...
            'MaximumNumPoints', max(numel(DataToAnalyze.BAT(1).AcousticSig_All),1) ); % Jamming signals
        hClutterTimePlot{nn}        = animatedline(h2, 'Color', colStruct.clutter,  'LineWidth', 2, ...
            'MaximumNumPoints', max(numel(DataToAnalyze.BAT(1).AcousticSig_All),1) ); % Clutter
        hConspsTimePlot{nn}         = animatedline(h2, 'Color', colStruct.conspsEchoes,  'LineWidth', 2, ...
            'MaximumNumPoints', max(numel(DataToAnalyze.BAT(1).AcousticSig_All),1) ); % Consps
    else % acoustics
        hCallsTimePlot{nn}          = animatedline(h2, 'Color', colStruct.ownCalls,  'LineWidth', 2, ...
            'MaximumNumPoints', max(size(DataToAnalyze.BAT(1).BatSonarEchosMat,2),1) ); % calls
        hWantedEchoesTimePlot{nn}   = animatedline(h2, 'Color', colStruct.wantedSignals,  'LineWidth', 2, ...
            'MaximumNumPoints', max(size(DataToAnalyze.BAT(1).BatSonarEchosMat,2),1) ); % wanted echoes echoes
        hMaskingSignalsTimePlot{nn} = animatedline(h2, 'Color',  colStruct.conspsCalls,  'LineWidth', 2, ...
            'MaximumNumPoints', max(size(DataToAnalyze.BAT(1).BatSonarEchosMat,2),1) ); % Jamming signals
        hClutterTimePlot{nn}        = animatedline(h2, 'Color', colStruct.clutter,  'LineWidth', 2, ...
            'MaximumNumPoints', max(size(DataToAnalyze.BAT(1).BatSonarEchosMat,2),1) ); % Clutter
        hConspsTimePlot{nn}         = animatedline(h2, 'Color', colStruct.conspsEchoes,  'LineWidth', 2, ...
            'MaximumNumPoints', max(size(DataToAnalyze.BAT(1).BatSonarEchosMat,2),1) ); % Consps
    end % if SimParams.AcousticsCalcFlag
    hCatchPreyTimes{nn} = animatedline(h2, 'Color', colStruct.catchPrey, 'Marker','+',...
            'MarkerSize', 8, 'linewidth', 1.2, 'LineStyle','none');
    hJammedTimes{nn}    = animatedline(h2, 'color', colStruct.masking, 'Marker','d',...
            'MarkerSize', 8, 'linewidth', 1.2, 'LineStyle','none');
    hCrushesTimes{nn}   = animatedline(h2, 'color', colStruct.crushes, 'Marker','x',...
            'MarkerSize', 10, 'linewidth', 1.8, 'LineStyle','none');
    hDetTimes{nn}       = animatedline(h2, 'color', colStruct.detections, 'Marker','*',...
            'MarkerSize', 8, 'linewidth', 1.5, 'LineStyle','none');
end

% text time
timeTxt   = text(h1, 0.82*h1.XLim(2), 1.05*h1.YLim(2), '0');
nPulseTxt = text(h1, 0.62*h1.XLim(2), 1.05*h1.YLim(2), '0');
grid(h2, 'minor')
h2.YTick = [0:20:140];

% %%% the Animation
% axes(h1);    
% AnimationTimeRate =  round(length(TimeVec)/(200/NumberOfBatsToPlot/max(1,NumberOfPreysToPlot) ));




%% THe Finding Esimated Postions Data
for nn = 1:NumberOfBats
    if DataToAnalyze.AllParams.SimParams.DetectObs
        obsFindsAll{nn} = obstacleDetectionSummary(DataToAnalyze.BAT(nn), DataToAnalyze.AllParams, 'Obs');
        obsFindsBatsindex(nn) = nn;
    end
    if DataToAnalyze.AllParams.SimParams.DetectPrey
        preyFindsAll{nn} = obstacleDetectionSummary(DataToAnalyze.BAT(nn), DataToAnalyze.AllParams, 'Prey');
        preyFindsBatsindex(nn) = nn;
    end
    if DataToAnalyze.AllParams.SimParams.DetectConsps
        conspsFindsAll{nn} = obstacleDetectionSummary(DataToAnalyze.BAT(nn), DataToAnalyze.AllParams, 'Cosnps');
        conspsFindsBatsindex(nn) = nn;
    end
end % for nn = 1:NumberOfBatsToPlot

%% ACoustics or not
batSigLvls(NumberOfBats) = struct('Calls',[], 'PreyEchoes',[], 'ConspsCalls',[], 'Clutter',[],  'ConspsEchoes',[]);
for kk = 1:NumberOfBats
    BATx       = DataToAnalyze.BAT(kk);    
    if SimParams.AcousticsCalcFlag && numel(BATx.AcousticSig_Calls) > 1
        batSigLvls(kk).Calls        = 20*log10( abs(BATx.AcousticSig_Calls)+1 ); 
        batSigLvls(kk).PreyEchoes   = 20*log10( abs(BATx.AcousticSig_PreyEchoes)+1 ); 
        batSigLvls(kk).ConspsCalls  = 20*log10( abs(BATx.AcousticSig_ConspsCalls)+1 ); 
        batSigLvls(kk).Clutter      = 20*log10( abs(BATx.AcousticSig_Clutter)+1 ); 
        batSigLvls(kk).ConspsEchoes = 20*log10( abs(BATx.AcousticSig_CospsEchoes)+1 ); 
    else
        batSigLvls(kk).Calls        = 10*log10( abs(BATx.BatSonarEchosMat(1,:))+1 ); 
        batSigLvls(kk).PreyEchoes   = 10*log10( abs(BATx.PreyEchosVec)+1 ); 
        batSigLvls(kk).ConspsCalls  = 10*log10( abs(BATx.AllInterPulses)+1 ); 
        batSigLvls(kk).Clutter      = 10*log10( abs(BATx.ObsEchosVec)+1) ; 
        batSigLvls(kk).ConspsEchoes = 10*log10( abs(BATx.Consps_EchosVec)+1 ); 
    end % Acoustics
end % for kbat

%% Starting Position of the wanted Bat
BATselected  = DataToAnalyze.BAT(DataToPlot.BatNumToPlot);
xSelectedBat = BATselected.BatX0;
ySelectedBat = BATselected.BatY0;

%% The animation
% plot(0,0,'Parent', axes2);
drawnow

% for kk = 1: AnimationTimeRate: length(TimeVec)-AnimationTimeRate+1
for kk = startXyIx: AnimationTimeRate: endXyIx-AnimationTimeRate+1

    for kBat = 1:NumberOfBatsToPlot
        CurMin     = kk;
        CurMax     = kk+ AnimationTimeRate -1;
        CurrentVec = CurMin:CurMax;
        BatToPlot  = BatNumToPlotVec(kBat);
        BATx       = DataToAnalyze.BAT(BatToPlot);
        
        % Acoustics
        if SimParams.AcousticsCalcFlag
            CurrentnTimeVec = int32( ((CurMin-1)*SampleTime*FsAcoustic + 1) : (CurMax*SampleTime*FsAcoustic) );
            ts = double(CurrentnTimeVec) * 1/FsAcoustic;
        else % Not-acoustics
            CurrentnTimeVec = CurMin:CurMax;
            ts = TimeVec(CurrentnTimeVec);
        end
        %% Clear Previous lines
        if AnimmateCurrentOnly
            clearpoints(hBatPos{kBat});
%             clearpoints(hSelectedBatPos);
            clearpoints(hCatchPreyPos{kBat})
            clearpoints(hBuzzStage{kBat})
            clearpoints(hObsManStage{kBat})
            clearpoints(hAVoitBatManStage{kBat})
            clearpoints(hCaveExitStage{kBat})
            clearpoints(hApprocahStage{kBat})
            clearpoints(hJammedPos{kBat})
            clearpoints(hCrushesPos{kBat})
            clearpoints(hClutterEstimates{kBat})
            clearpoints(hClutterFinds{kBat})
            clearpoints(hPreyEstimates{kBat});
            clearpoints(hPreyFinds{kBat});
%             clearpoints(hConspsFinds{kBat});
        end % if AnimmateCurrentOnly

        %%% Position of the Bats %%%
        xBatPos = BATx.xBatPos(CurrentVec);
        yBatPos = BATx.yBatPos(CurrentVec);
        %         set(0, 'CurrentFigure', f1);
        
        addpoints( hBatPos{kBat}, xBatPos, yBatPos );
%         addpoints( hSelectedBatPos, xSelectedBat, ySelectedBat );
        
        if BatToPlot ==  DataToPlot.BatNumToPlot 
           clearpoints(hSelectedBatPos);
           xSelectedBat = xBatPos;
           ySelectedBat = yBatPos;
           addpoints( hSelectedBatPos, xSelectedBat, ySelectedBat ); 
           drawnow;
        end
        %%% Catches Times %%%
        CatchesTimes = BATx.CatchPreyTimes./SampleTime;
        CatchesInd = find( (CatchesTimes >= CurMin) & (CatchesTimes <= CurMax));

        if ~isempty(CatchesInd)
            CatchPreyPos = BATx.CatchPreyPos(CatchesInd,:);
            xCactches = CatchPreyPos(:,1);
            yCactches = CatchPreyPos(:,2);
            addpoints( hCatchPreyPos{kBat}, xCactches, yCactches);
        end% if ~isempty(find)

        %%% STAGES Plot By color
        %%

        FullFlightStages = BATx.InterReportStrctOnLine.FullFlightStages(CurMin:CurMax);

        %%% Stage "Buzz" %%%
        BuzzIndex1 = strfind(FullFlightStages,'Buzz');
        BuzzIndex = find(~(cellfun('isempty', BuzzIndex1)));
        
        if ~isempty(ColorObsMan)
            addpoints(hBuzzStage{kBat}, xBatPos(BuzzIndex),yBatPos(BuzzIndex)) % the position when Approach Manuvering
        end % if ~isempty(HuntingIndex)

        %%% Stage "ObstacleManuver" %%%
        ObsManIndex1 = strfind(FullFlightStages,'ObstacleManuver');
        ObsManIndex = find(~(cellfun('isempty', ObsManIndex1)));
        if ~isempty(ObsManIndex1)
            addpoints(hObsManStage{kBat}, xBatPos(ObsManIndex),yBatPos(ObsManIndex)) % the position when Approach Manuvering
        end % if ~isempty(HuntingIndex)

        %%% Stage "Approach" %%%
        ApproachIndex1 = strfind(FullFlightStages,'Approach');
        ApproachIndex = find(~(cellfun('isempty', ApproachIndex1)));
        if ~isempty(ApproachIndex)
            addpoints(hApprocahStage{kBat}, xBatPos(ApproachIndex),yBatPos(ApproachIndex)) % the position when Approach Manuvering
        end % if ~isempty(HuntingIndex)
         
        %%% Stage AvoidBat %%%%
        
        AvoidBatIndex = find(strcmp(FullFlightStages,'AvoidBatMan'));
        if ~isempty(AvoidBatIndex)
             addpoints(hAVoitBatManStage{kBat}, xBatPos(AvoidBatIndex),yBatPos(AvoidBatIndex)) %  the position whenManuvering
        end % if ~isempty(HuntingIndex)


         %%% Stage "caveExit" %%%
        exitIndex1 = strfind(FullFlightStages,'CaveExit');
        exitIndex = find(~(cellfun('isempty', exitIndex1)));
        if ~isempty(exitIndex)
            addpoints(hCaveExitStage{kBat}, xBatPos(exitIndex),yBatPos(exitIndex)) % the position when Approach Manuvering
        end % if ~isempty(HuntingIndex)
        

        %% The Findings of Prey in each sequence
        if DataToAnalyze.AllParams.SimParams.DetectPrey
            batIx = preyFindsBatsindex == BatToPlot;
            preyFindsIx = find( ([preyFindsAll{batIx}.nTime] >= CurMin) & ([preyFindsAll{batIx}.nTime] <= CurMax));
            if ~isempty(preyFindsIx)
                xPreyReal = [preyFindsAll{batIx}(preyFindsIx).xObsReal];
                yPreyReal = [preyFindsAll{batIx}(preyFindsIx).yObsReal];
                addpoints(hPreyFinds{kBat}, xPreyReal, yPreyReal);

                xPreyEstimated = [preyFindsAll{batIx}(preyFindsIx).xObsEstimated];
                yPreyEstimated = [preyFindsAll{batIx}(preyFindsIx).yObsEstimated];
                %             addpoints(hPreyEstimates{kBat}, xPreyEstimated, yPreyEstimated);
            end % if ~isempty(clutterIx)
        end   %     if DataToAnalyze.AllParams.SimParams.DetectPrey



        %%
        %%% JammedTimes %%%
%         JammedTimes = ...
%             BATx.InterReportStrctOnLine.TotalInterferenceTimes / SampleTime;

        switch DetectedObj
            case 'Prey'
                JammedTimes = BATx.InterReportStrctOnLine.TotalInterferenceTimes / SampleTime;
                DetTimes    = [BATx.PreyFindsStruct.DetectedTimes];
                wantedLvls  = batSigLvls(kBat).PreyEchoes;
            case 'Obs'
                JammedTimes = BATx.InterReportStrctOnLine.obsTotalInterferenceTimes / SampleTime;
                DetTimes    = [BATx.ObsFindsStruct.DetectedTimes];
                wantedLvls  = batSigLvls(kBat).Clutter;
            case 'Conspecific'
                JammedTimes = BATx.InterReportStrctOnLine.conspsTotalInterferenceTimes / SampleTime;
                DetTimes    = [BATx.Consps_FindsStruct.DetectedTimes];
                wantedLvls  = batSigLvls(kBat).ConspsEchoes;
        end % switch

        JammsInd      = find( (JammedTimes >= CurMin) & (JammedTimes <= CurMax));
        DetectionsInd = find( (DetTimes >= CurMin) & (DetTimes <= CurMax));

%         if ~isempty(JammsInd)
%             xJammedPos = BATx.xBatPos(round(JammedTimes(JammsInd)));
%             yJammedPos = BATx.yBatPos(round(JammedTimes(JammsInd)));
%             addpoints(hJammedPos{kBat}, xJammedPos, yJammedPos);
% 
%         end % if ~isempty(JammedTimes)
        
        %% Crushes
        %%% obs
        crushIx = find( (BATx.CrushesObsnTimes >= CurMin) & (BATx.CrushesObsnTimes <= CurMax));
        if ~isempty(crushIx)
%         if BATx.CrushesObsNum > 0
            xCrushPos = BATx.xBatPos(round(BATx.CrushesObsnTimes(crushIx)));
            yCrushPos = BATx.yBatPos(round(BATx.CrushesObsnTimes(crushIx)));
            addpoints(hCrushesPos{kBat}, xCrushPos, yCrushPos);

        end % if ~isempty(crushIx)
        
        %%% consp
        crushIx = find( (BATx.CrushesConspsnTimes >= CurMin) & (BATx.CrushesConspsnTimes <= CurMax));
        if ~isempty(crushIx)
%         if BATx.CrushesObsNum > 0
            xCrushPos = BATx.xBatPos(round(BATx.CrushesConspsnTimes(crushIx)));
            yCrushPos = BATx.yBatPos(round(BATx.CrushesConspsnTimes(crushIx)));
            addpoints(hCrushesPos{kBat}, xCrushPos, yCrushPos);

        end % if ~isempty(crushIx)



        %% update time of simulation
        timeTxt.String = ['time(sec): ', num2str(ts(1),3)];
        %find the first pulse number in ix
        if BatToPlot ==  DataToPlot.BatNumToPlot
            nPulse = find([BATx.TransmittedPulsesStruct.StartPulseTime] >= CurMin,1);
            nPulseTxt.String = ['PulseNum: ', num2str(nPulse)]; 
        end % if BatToPlot ==  DataToPlot.BatNumToPlot

        %% The Clutter  -Only for the Selected Bat
        if BatToPlot ==  DataToPlot.BatNumToPlot
            if DataToAnalyze.AllParams.SimParams.DetectObs
                batIx = obsFindsBatsindex == BatToPlot;
                clutterIx = find( ([obsFindsAll{batIx}.nTime] >= CurMin) & ([obsFindsAll{batIx}.nTime] <= CurMax));
                if ~isempty(clutterIx)
                    xObsReal = [obsFindsAll{batIx}(clutterIx).xObsReal];
                    yObsReal = [obsFindsAll{batIx}(clutterIx).yObsReal];
                    % addpoints(hClutterFinds{kBat}, xObsReal, yObsReal);

                    xObsEstimated = [obsFindsAll{batIx}(clutterIx).xObsEstimated];
                    yObsEstimated = [obsFindsAll{batIx}(clutterIx).yObsEstimated];
                    addpoints(hClutterEstimates{kBat}, xObsEstimated, yObsEstimated);
                end % if ~isempty(clutterIx)
            end %if DataToAnalyze.AllParams.SimParams.DetectObs
        end %if BatToPlot ==  DataToPlot.BatNumToPlot

        %% The Findings of Consps in each sequence - Only for the Selected Bat
        if BatToPlot ==  DataToPlot.BatNumToPlot
            if DataToAnalyze.AllParams.SimParams.DetectConsps

                clearpoints(hConspsFinds{kBat});
                batIx = conspsFindsBatsindex == BatToPlot;
                conspsFindsIx = find( ([conspsFindsAll{batIx}.nTime] >= CurMin) & ([conspsFindsAll{batIx}.nTime] <= CurMax));
                if ~isempty(conspsFindsIx)
                    xConspsReal = [conspsFindsAll{batIx}(conspsFindsIx).xObsReal];
                    yConspsReal = [conspsFindsAll{batIx}(conspsFindsIx).yObsReal];
                    %                 addpoints(hConspsFinds{kBat}, xConspsReal, yConspsReal);
                    %                 %% Tobe Fixed

                    xConspsEstimated = [conspsFindsAll{batIx}(conspsFindsIx).xObsEstimated];
                    yConspsEstimated = [conspsFindsAll{batIx}(conspsFindsIx).yObsEstimated];
                    addpoints(hConspsFinds{kBat}, xConspsEstimated, yConspsEstimated);
                    %             addpoints(hPreyEstimates{kBat}, xPreyEstimated, yPreyEstimated);
                end % if ~isempty(clutterIx)
            end %         if DataToAnalyze.AllParams.SimParams.DetectCons
             
            %%%%% The jammed objects of the selected bat
            if ~isempty(JammsInd)
                
                xJammedPos = BATx.xBatPos(round(JammedTimes(JammsInd)));
                yJammedPos = BATx.yBatPos(round(JammedTimes(JammsInd)));
                addpoints(hJammedPos{kBat}, xJammedPos, yJammedPos);

            end % if ~isempty(JammsInd)

        end % if BatToPlot

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Time Plot
        %         hCallsTimePlot{NumberOfBatsToPlot} = animatedline(h2, 'Color', 'k'); % calls
        %         hWantedEchoesTimePlot{NumberOfBatsToPlot} = animatedline(h2, 'Color', 'g'); % wanted echoes echoes
        %         hMaskingSignalsTimePlot{NumberOfBatsToPlot} = animatedline(h2, 'Color', 'r'); % Jamming signals
        
        % plot time plot only for the BatToPlot
        if BatToPlot ==  DataToPlot.BatNumToPlot
%             clearpoints(hSelectedBatPos);
            try
%                 %%% ACoustics or not
%                 if SimParams.AcousticsCalcFlag && numel(BATx.AcousticSig_Calls) > 1
%                     addpoints(hCallsTimePlot{kBat},          ts, 20*log10( abs(BATx.AcousticSig_Calls(CurrentnTimeVec) +1)) ) % the position when Approach Manuvering
%                     addpoints(hWantedEchoesTimePlot{kBat},   ts, 20*log10( abs(BATx.AcousticSig_PreyEchoes(CurrentnTimeVec) +1)) )
%                     addpoints(hMaskingSignalsTimePlot{kBat}, ts, 20*log10( abs(BATx.AcousticSig_ConspsCalls(CurrentnTimeVec) +1)) )
%                     addpoints(hClutterTimePlot{kBat},        ts, 20*log10( abs(BATx.AcousticSig_Clutter(CurrentnTimeVec) +1)) )
%                     addpoints(hConspsTimePlot{kBat},         ts, 20*log10( abs(BATx.AcousticSig_CospsEchoes(CurrentnTimeVec) +1)) )
%                 else
%                     addpoints(hCallsTimePlot{kBat},          tt, 10*log10( abs(BATx.BatSonarEchosMat(1,currIx) +1)) ) % the position when Approach Manuvering
%                     addpoints(hWantedEchoesTimePlot{kBat},   tt, 10*log10( abs(BATx.PreyEchosVec(currIx) +1)) )
%                     addpoints(hMaskingSignalsTimePlot{kBat}, tt, 10*log10( abs(BATx.AllInterPulses(currIx) +1)) )
%                     addpoints(hClutterTimePlot{kBat},        tt, 10*log10( abs(BATx.ObsEchosVec(currIx) +1)) ) 
%                     addpoints(hConspsTimePlot{kBat},         tt, 10*log10( abs(BATx.Consps_EchosVec(currIx) +1)) ) 
%                 end % Acoustics
                addpoints(hCallsTimePlot{kBat},          ts, batSigLvls(BatToPlot).Calls(CurrentnTimeVec) ) % the position when Approach Manuvering
                addpoints(hMaskingSignalsTimePlot{kBat}, ts, batSigLvls(BatToPlot).ConspsCalls(CurrentnTimeVec) )
                addpoints(hConspsTimePlot{kBat},         ts, batSigLvls(BatToPlot).ConspsEchoes(CurrentnTimeVec) )
                addpoints(hClutterTimePlot{kBat},        ts, batSigLvls(BatToPlot).Clutter(CurrentnTimeVec) )
                addpoints(hWantedEchoesTimePlot{kBat},   ts, batSigLvls(BatToPlot).PreyEchoes(CurrentnTimeVec) )
                
                %%% Moving Detection by Time
%                 xlim(h2,[CurMin*SampleTime-1, CurMax*SampleTime])

                %%% Capture Times
                if ~isempty(CatchesInd)
                    currCatchTimes = CatchesTimes(CatchesInd) * SampleTime;
                    addpoints( hCatchPreyTimes{kBat}, currCatchTimes, 80*ones(size(currCatchTimes)) )
                end% if ~isempty(find)

                %%% jAMMED tIMES
                if ~isempty(JammsInd)
                    currJammedTimes = JammedTimes(JammsInd) * SampleTime;

                    addpoints(hJammedTimes{kBat}, currJammedTimes, 80*ones(size(JammsInd)) )
                end % if ~isempty(JammedTimes)
                
                %%% Detection Times - dont plot for clutter
                if ~isempty(DetectionsInd) && ~strcmp(DetectedObj, 'Obs') %'Prey'; % 'Prey'; 'Obs'   
                    
                    if SimParams.AcousticsCalcFlag
                        ixx = round(DetTimes(DetectionsInd)*FsAcoustic*SampleTime);
                    else
                        ixx = round(DetTimes(DetectionsInd));
                    end
                    % Original
                    currDetTimes = DetTimes(DetectionsInd) * SampleTime;
                    currDetLvl   =  wantedLvls(ixx); % repmat(currDetLvl, 1, numel(currDetTimes)) 
                    % addpoints(hDetTimes{kBat}, currDetTimes, currDetLvl)

                    % after Fade - tale the signal Level after rise time
                    % (some degree)
                    scanSamples = 500;
                    mTimes = currDetTimes;
                    mLvl   = currDetLvl;
                    for iPk = 1:numel(ixx)
                        [mLvl(iPk), mTimes(iPk)] = max(wantedLvls(ixx(iPk):ixx(iPk)+scanSamples));
                    end
                    currDetTimes = currDetTimes + mTimes/FsAcoustic;
                    currDetLvl = mLvl;
                    addpoints(hDetTimes{kBat}, currDetTimes, currDetLvl)
                end % if ~isempty(JammedTimes)


                %%% Crushes
                % obs
                crushIx = find( (BATx.CrushesObsnTimes >= CurMin) & (BATx.CrushesObsnTimes <= CurMax));
                if ~isempty(crushIx)
                    %         if BATx.CrushesObsNum > 0
                    addpoints(hCrushesTimes{kBat}, BATx.CrushesObsnTimes(crushIx)*SampleTime, 100*ones(size(crushIx)) );
                end % if ~isempty(crushIx)
                % consps
                crushIx = find( (BATx.CrushesConspsnTimes >= CurMin) & (BATx.CrushesConspsnTimes <= CurMax));
                if ~isempty(crushIx)
                    %         if BATx.CrushesObsNum > 0
                    addpoints(hCrushesTimes{kBat}, BATx.CrushesConspsnTimes(crushIx)*SampleTime, 100*ones(size(crushIx)) );
                end % if ~isempty(crushIx)
            catch
                warning(['TimeVec Error', num2str(ts(1))] )
            end % try
        end % if kBat ==  DataToPlot.BatNumToPlot


        drawnow
         if SaveVideoFlag
            CurrFrame = getframe(gcf);
            writeVideo(VideoObj,CurrFrame)
        end % if SaveVideoFlag

    end % for kBat = 1:NumberOfBatsToPlot

    for kPrey = 1:NumberOfPreysToPlot

        %         clearpoints( hPreyPos{kPrey});
        %         drawnow
        
        PreyToPlot = PreyNumToPlotVec(kPrey);
        xPreyPos = DataToAnalyze.PREY(PreyToPlot).xPrey;
        yPreyPos = DataToAnalyze.PREY(PreyToPlot).yPrey;
        %         set(0, 'CurrentFigure', f1);
        if AnimmateCurrentOnly
            clearpoints(hPreyPos{kPrey});
            clearpoints(hPreyJamPos{kPrey});
        end %  if AnimmateCurrentOnly
        
        addpoints( hPreyPos{kPrey}, xPreyPos(kk:kk+ AnimationTimeRate -1), yPreyPos(kk:kk+ AnimationTimeRate -1));

        drawnow
        if SaveVideoFlag
            CurrFrame = getframe(gcf);
            writeVideo(VideoObj,CurrFrame)
        end % if SaveVideoFlag

        %%%% THe jamming Tx times of Preys
        if DataToPlot.JamPreyPosFlag

            nJamTx = round(DataToAnalyze.PREY(PreyToPlot).JammingnTimes);
            ixJam = find( (nJamTx >= CurMin) & (nJamTx <= CurMax));
            
            addpoints( hPreyJamPos{kPrey}, DataToAnalyze.PREY(PreyToPlot).xPrey(nJamTx(ixJam)), DataToAnalyze.PREY(PreyToPlot).yPrey(nJamTx(ixJam)));
%             plot(h1, DataToAnalyze.PREY(PreyToPlot).xPrey(nJamTx(ixJam)), DataToAnalyze.PREY(PreyToPlot).yPrey(nJamTx(ixJam)), 'Color', 'r', 'Marker', 'x', ...
%                 'MarkerSize',9, 'LineStyle', 'none') % the position of the Prey
            

        end % if DataToPlot.JamPreyPosFlag

    end % for kPrey = 1:NumberOfBatsToPlot

end % for kk = 1:length(xBatPos)

if SaveVideoFlag
    close(VideoObj);
end % if SaveVideoFlag





end % function