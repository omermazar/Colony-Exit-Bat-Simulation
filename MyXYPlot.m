function [] = MyXYPlot(HxyPlot,DataToAnalyze,DataToPlot,...
    SimParams, TerrainParams, Terrain, ZoomFlag, ZoomIDn)
% Changed NOV2021 - New Terrain
% Cahnge 17May2025 - FullFlightStages is by ManuverType and not Manuver stage

%% INIT
h1 = HxyPlot;
% hold(h1,'on');
NumberOfBats = DataToAnalyze.AllParams.SimParams.TotalBatsNumber;
NumberOfPreys= DataToAnalyze.AllParams.SimParams.TotalPreysNumber;

SampleTime =  DataToAnalyze.AllParams.SimParams.SampleTime;
MaxTime =  DataToAnalyze.AllParams.SimParams.SimulationTime;
TimeVec = [0:SampleTime:MaxTime]; % Timee Vector
xyResolution = DataToAnalyze.AllParams.SimParams.xyResolution;



%%% Which Bats to plot %%%
if DataToPlot.PlotAllBatsFlag 
    NumberOfBatsToPlot = DataToAnalyze.AllParams.SimParams.TotalBatsNumber;
    BatNumToPlotVec = 1:NumberOfBatsToPlot;
else %if DataToPlot.PlotAllBatsFlag 
    NumberOfBatsToPlot =1;
    BatNumToPlotVec = [DataToPlot.BatNumToPlot];  
end % if DataToPlot.PlotAllBatsFlag 

%%% Which Preys to plot %%%
if DataToPlot.PlotAllPreysFlag 
    NumberOfPreysToPlot = DataToAnalyze.AllParams.SimParams.TotalPreysNumber;
    PreyNumToPlotVec = 1:NumberOfPreysToPlot;
else %if DataToPlot.PlotAllBatsFlag 
    NumberOfPreysToPlot =1;
    PreyNumToPlotVec = [DataToPlot.PreyNumToPlot];  
end % if Data

if NumberOfPreys >0
    PREY(NumberOfPreys) = struct('xPrey',[], 'yPrey',[]);
end % if NumberOfPreys >0

StartTimeIndex = 1; % 1;
EndTimeIndex = MaxTime/SampleTime+1;

%%% Zoom %%%
% if DataToPlot.TimeFilterFlag && ~ZoomFlag %
if DataToPlot.TimeFilterFlag
    cla(h1);
    StartTimeIndex = max(DataToPlot.StartTimeFilter./SampleTime , 1);
    EndTimeIndex = DataToPlot.EndTimeFilter./SampleTime;
end % if 

IDn = StartTimeIndex:EndTimeIndex;

if ZoomFlag
    cla(h1);
    IDn = ZoomIDn;
end % ZoomFlag

% StartTimeIndex
TimeVec = TimeVec(IDn);

% for nPrey = 1:NumberOfPreys
% %         PREY(nPrey).xPrey = DataToAnalyze.PREY(nPrey).xPrey(IDn);
% %         PREY(nPrey).yPrey = DataToAnalyze.PREY(nPrey).yPrey(IDn);
%         
% end % for n = 1:NumberOfPreys

PREY = DataToAnalyze.PREY;


%%% Setting Colors
TotalBatsNumber = DataToAnalyze.AllParams.SimParams.TotalBatsNumber;
TotalPreysNumber = DataToAnalyze.AllParams.SimParams.TotalPreysNumber;

[colStruct] = myColorScheme(SimParams.TestMode, TotalBatsNumber, TotalPreysNumber);

% kColor1 = linspace(0.5,1, TotalBatsNumber);
% ColorOfBatPos = zeros(TotalBatsNumber,3);
% % ColorOfBatPos(:,3) = kColor1; % Blue
% ColorOfBatPos = lines(NumberOfBatsToPlot); % the colors of the lines for bats

% kColor2 = linspace(0.5,1, TotalPreysNumber);
% ColorOfPreysPos = zeros(TotalPreysNumber,3);
% ColorOfPreysPos(:,3) = kColor2;
% ColorOfPreysPos(:,1) = kColor2; % magenta

ColorForagingMan = [ 0    0.4470    0.7410]; % 'b'; colStruct.ForagingMan;
ColorHuntingMan  = colStruct.Approach; % [ 0 1 0];  % Green
ColorBuzzMan     = colStruct.BuzzMan; %[ 0 0 0 ]; % black
ColorObsMan      = 'y'; colStruct.ObsMan; %[0 0.7 0]; % Green
ColorFindingPrey = colStruct.preyFinding;% [1 0 1]; % magenta
ColorAvoidBatMan = 'y'; colStruct.AvoidBatMan; %'y';
ColorFollow      = 'g';

%%% THE PLOTS %%%
%%
% Plot the terrain
if DataToPlot.TerrainFlag
    %Plot the terrain and the objects inside it
    sz = size(Terrain,1);
%     sz=sz(1);
    for i=1:+2:sz-1
        plot(h1, Terrain(i,:)*xyResolution,Terrain(i+1,:)*xyResolution,'b-', 'DisplayName', 'borders')
    end
end % HDataToPlot.TerrainFlag
% % % if DataToPlot.TerrainFlag
% % %     [Xcor, Ycor] = meshgrid(0:SimParams.xyResolution:TerrainParams.Xmax,...
% % %         0:SimParams.xyResolution:TerrainParams.Ymax);
% % %     if ~strcmp(DataToAnalyze.AllParams.TerrainParams.Type, 'Free Space')
% % %         contour(h1,Xcor,Ycor,Terrain);
% % %     end % if
% % % end % HDataToPlot.TerrainFlag

%% Bats PLOTS:  Position, SOnar Places, Manuvers Etc
for kBat = 1:NumberOfBatsToPlot
    BatToPlot = BatNumToPlotVec(kBat);
    xBatPos = DataToAnalyze.BAT(BatToPlot).xBatPos(IDn);
    yBatPos = DataToAnalyze.BAT(BatToPlot).yBatPos(IDn);
    BatX0 = DataToAnalyze.BAT(BatToPlot).BatX0;
    BatY0 = DataToAnalyze.BAT(BatToPlot).BatY0;
   
    % BatPosition Plot
    if DataToPlot.BatPosFlag
        %%% The Focal Bat
        ii = min(BatToPlot, size(colStruct.BatPos,1));
        if kBat == DataToPlot.BatNumToPlot & TotalBatsNumber > 1
            plot(h1,xBatPos,yBatPos,'color', colStruct.SelectedBatPos ,'linewidth',4, 'DisplayName', ['FocalBat', num2str(kBat)])
            plot(h1,xBatPos,yBatPos,'Color', colStruct.BatPos(BatToPlot,:) ,'MarkerFaceColor', colStruct.BatPos(BatToPlot,:), ... 
                'linewidth',2, 'DisplayName', ['Bat', num2str(kBat)])
        else
            ii = min(BatToPlot, size(colStruct.BatPos,1));
            plot(h1,xBatPos,yBatPos, 'color',colStruct.BatPos(ii,:) ,'linewidth',2, 'DisplayName', ['Bat', num2str(kBat)])
        end % if kBat == DataToPlot.BatNumToPlot
        plot(h1,BatX0, BatY0,'S','color', colStruct.BatPos(ii,:), 'MarkerSize',8,'linewidth',2, 'DisplayName', ['StartBat', num2str(kBat)]); 
        
        %%% Add Crushes 
%         plot(h1, xBatPos(DataToAnalyze.BAT(BatToPlot).CrushesObsnTimes), yBatPos(DataToAnalyze.BAT(BatToPlot).CrushesObsnTimes), ...
%             'xr','MarkerSize', 12, 'linewidth', 4, 'DisplayName', 'ObsCrush')
        
    end % HDataToPlot.BatPosFlag
    
    % Bat SonarPlaces Plot
    if DataToPlot.SonarPosFlag
        SonarTimesVec = DataToAnalyze.BAT(BatToPlot).SonarPulses(IDn);
        FindsVec = find(SonarTimesVec);
        plot(h1,xBatPos(FindsVec),yBatPos(FindsVec),'*-', 'color', colStruct.EcholocatePos,  'DisplayName', 'SonarPos') % the position when sonar
%         plot(h1,xBatPos(find(SonarTimesVec)),yBatPos(find(SonarTimesVec)),'g-') % the position when sonar
    end % HDataToPlot.SonarPosFlag
    
    % Bat Manuvers to avoid Obstacles, to Hunt and the catches positions of the prey
    if DataToPlot.ManPosFlag
        %%%% new 17May2025
        nPulses = DataToAnalyze.BAT(BatToPlot).NumOfTimesSonarTransmits;
        typeCell = {DataToAnalyze.BAT(BatToPlot).ManueverCmdStruct(1:nPulses).ManueverType};
        PulseTimesVec = [ DataToAnalyze.BAT(BatToPlot).TransmittedPulsesStruct.StartPulseTime];
        FullFlightStages = Pulses2Flight(typeCell, PulseTimesVec, nPulses);
        % clean arrray
        ix = find(cellfun(@isnumeric, FullFlightStages));
        [FullFlightStages{ix}] = deal('');

        %%% old
            % FullFlightStages =  DataToAnalyze.BAT(BatToPlot).InterReportStrctOnLine.FullFlightStages;
        %%%%
        uStages =  unique(FullFlightStages(~cellfun('isempty',FullFlightStages))); %unique(FullFlightStages);
        % remove empty elments in the cell
        uStages =uStages(~cellfun('isempty',uStages));
        
%         HuntingIndexCell = strfind(DataToAnalyze.BAT(BatToPlot).ManueverTypeCell,'Hunting');
        ApproachIndexCell = strfind(FullFlightStages,'Approach');
        ApproachIndex = find(~(cellfun('isempty', ApproachIndexCell)));
        if ~isempty(ApproachIndex)
            plot(h1,xBatPos(ApproachIndex),yBatPos(ApproachIndex), '.', 'color', ColorHuntingMan ,'linewidth', 12, 'DisplayName', 'Approach') % the position whenManuvering
        end % if ~isempty(HuntingIndex)
        
        % avoid obstacle
        ObsManIndexCell = strfind(FullFlightStages,'ObsMan'); % strfind(FullFlightStages,'ObstacleManuver');
        ObsManIndex = find(~(cellfun('isempty', ObsManIndexCell)));
        if ~isempty(ObsManIndex)
            plot(h1,xBatPos(ObsManIndex),yBatPos(ObsManIndex),'.', 'color', ColorObsMan ,'MarkerSize', 9, 'DisplayName', 'ObstacleMan') % the position whenManuvering
        end % if ~isempty(HuntingIndex)
        
        BuzzIndexCell = strfind(FullFlightStages,'Buzz');
        BuzzIndex = find(~(cellfun('isempty', BuzzIndexCell)));
        if ~isempty( BuzzIndex)
            % plot(h1,xBatPos( BuzzIndex),yBatPos( BuzzIndex),'.', 'color', ColorBuzzMan ,'MarkerSize',7, 'DisplayName', 'Buzz') % the position whenManuvering
        end % if ~isempty(HuntingIndex)
        
        AvoidBatIndex = find(strcmp(FullFlightStages,'AvoidBatMan'));
        if ~isempty(AvoidBatIndex)
            plot(h1,xBatPos(AvoidBatIndex),yBatPos(AvoidBatIndex),'.', 'color', ColorAvoidBatMan ,'MarkerSize', 9, 'DisplayName', 'ConspAvoid') % the position whenManuvering
        end % if ~isempty(HuntingIndex)

        caveExitIndex = find(strcmp(FullFlightStages,'ExitMan')); % find(strcmp(FullFlightStages,'CaveExit'));
        if ~isempty(caveExitIndex)
            plot(h1,xBatPos(caveExitIndex),yBatPos(caveExitIndex),'.', 'color', colStruct.CaveExit ,'MarkerSize', 9, 'DisplayName', 'CaveExit') % the position whenManuvering
        end % if ~isempty(HuntingIndex)

        % Follow Wall
        followIndex = find(strcmp(FullFlightStages,'FollowWall'));
        if ~isempty(followIndex)
            plot(h1,xBatPos(followIndex),yBatPos(followIndex),'.', 'color', ColorFollow ,'MarkerSize', 9, 'DisplayName', 'CaveExit') % the position whenManuvering
        end % if ~isempty(HuntingIndex)
        
        % Foraging
        forageIndex = find(strcmp(FullFlightStages,'Foraging'));
        if ~isempty(followIndex)
            plot(h1,xBatPos(forageIndex),yBatPos(forageIndex),'.', 'color', ColorForagingMan ,'MarkerSize', 9, 'DisplayName', 'CaveExit') % the position whenManuvering
        end % if ~isempty(HuntingIndex)
    end % DataToPlot.ManPosFlag
    
  % The Places Where Jamming accured
    DataToPlot.JammingPlacesFlag = 1; 
    if DataToPlot.JammingPlacesFlag
        JammedTimes = ...
            DataToAnalyze.BAT(BatToPlot).InterReportStrctOnLine.TotalInterferenceTimes...
             / SampleTime ;
        if ~isempty(JammedTimes)
            JammedBatPos = [DataToAnalyze.BAT(BatToPlot).xBatPos(round(JammedTimes)) ; ...
                DataToAnalyze.BAT(BatToPlot).yBatPos(round(JammedTimes))];
            plot(h1, JammedBatPos(1,:), JammedBatPos(2,:),'db','MarkerSize',2,'color', colStruct.Masking,'linewidth',4, 'DisplayName', 'JammedEcho')
        end % if ~isempty(JammedTimes)
    end % if DataToPlot.JammingPlacesFlag
    
    
    % the Obstacles found by bat
    if DataToPlot.TerrainFindFlag && BatToPlot == DataToPlot.BatNumToPlot % plot only for the foacl bat
        xObsPos = DataToAnalyze.BAT(BatToPlot).xObsFinds;
        yObsPos = DataToAnalyze.BAT(BatToPlot).yObsFinds;
        try
            % plot(h1,xObsPos,yObsPos,'o', 'Color', colStruct.Terrain, ...
            %      'LineWidth',1, 'MarkerSize', 6) % the position of the obsitcles
        catch
            warning('Problem plotting Terrain finds!!! (line 187)')
        end
        % The Obstacles in the bat point of view
        estObsFlag = true;
        if estObsFlag 
            [obsFindsAll] = obstacleDetectionSummary(DataToAnalyze.BAT(BatToPlot), DataToAnalyze.AllParams, 'Obs');
            plot([obsFindsAll.xObsEstimated], [obsFindsAll.yObsEstimated], 's', 'MarkerFaceColor', colStruct.obsEstimate, ...
                'MarkerEdgeColor', colStruct.obsEstimate, 'MarkerSize', 4, 'DisplayName', 'ObstacleDetections')
        end % if estObsFlag 
    end % HDataToPlot.SonarPosFlag
    
    % the Prey Movement
    
    
    % The positions of the prey found by the bat
    if DataToPlot.PreyFindsFlag
        xPreyFinds = nonzeros(DataToAnalyze.BAT(BatToPlot).xPreyFinds) ;
        yPreyFinds = nonzeros(DataToAnalyze.BAT(BatToPlot).yPreyFinds);
        nTimesPreyFinds = nonzeros(DataToAnalyze.BAT(BatToPlot).nTimesPreyFinds);
        %%% dealing with zoom
        StartZoomIndex = min(IDn);
        StopZoomIndex = max(IDn);
        tt1= find((nTimesPreyFinds >= StartZoomIndex) & (nTimesPreyFinds<= StopZoomIndex));
        xPreyFinds = xPreyFinds(tt1);
        yPreyFinds = yPreyFinds(tt1);
        tPreyFinds = nTimesPreyFinds(tt1);
        %%% plot
        plot(h1,xPreyFinds,yPreyFinds,'+', 'color', ColorFindingPrey, 'DisplayName', 'PreyDetections') % the position of the Prey
%         plot(h1,DataToAnalyze.BAT(BatToPlot).xBatPos(tPreyFinds),DataToAnalyze.BAT(BatToPlot).yBatPos(tPreyFinds), 'Marker','+',...
%             'MarkerSize',6, 'color', ColorFindingPrey, 'LineStyle','none') % the position of th bat  when finding prey
        %         end % if NumberOfTimesDetectedPrey >0
    end % if DataToPlot.PreyFindsFlag
    
    % The Catches of the Prey positions
    if DataToPlot.CatchesPosFlag
        CatchPreyTimes = DataToAnalyze.BAT(BatToPlot).CatchPreyTimes./SampleTime;
        CatchPreyPos = DataToAnalyze.BAT(BatToPlot).CatchPreyPos;
        
        if ~isempty(CatchPreyPos)
            for nCatch = 1: DataToAnalyze.BAT(BatToPlot).NumberOfCatches
                Index = find(IDn == round(CatchPreyTimes(nCatch)),1); % Check if the catch in ZOOM
                if ~isempty(Index) 
                    plot(h1, CatchPreyPos(nCatch,1), CatchPreyPos(nCatch,2), '+', 'color', colStruct.catchPrey, 'MarkerSize',16,'linewidth',4, 'DisplayName', 'Capture')
                end % if ~isempty(Index)
            end % for nCatch
        end % if ~isempty(CatchPreyPos)
    end % DataToPlot.CatchesPosFlag

end %for kBat = 1:NumberOfBatsToPlot   

%% Prey Plots
if DataToPlot.PreyPositionFlag
%         Color = ['b','c','g', 'r', 'k', 'y', 'm'];
        for kPrey = 1:NumberOfPreysToPlot
            
            plot(h1,PREY(PreyNumToPlotVec(kPrey)).xPrey(IDn), PREY(PreyNumToPlotVec(kPrey)).yPrey(IDn), 'Color',colStruct.PreysPos(kPrey,:), 'Marker','.', ...
                'LineStyle','none', 'DisplayName', ['Prey', num2str(kPrey)] ) % the position of the Prey
            
            plot(h1,DataToAnalyze.PREY(PreyNumToPlotVec(kPrey)).PreyX0*xyResolution, ...
                DataToAnalyze.PREY(PreyNumToPlotVec(kPrey)).PreyY0*xyResolution,'S', ...
                'Color',colStruct.PreysPos(kPrey,:), 'MarkerSize',6,'linewidth',2, 'DisplayName', ['StartPrey', num2str(kPrey)] )
      
        end % for kPrey
end %if DataToPlot.PreyPositionFlag

%%%% THe jamming Tx times of Preys
 if DataToPlot.JamPreyPosFlag
        
     for kPrey = 1:NumberOfPreysToPlot
         ixJam = round(PREY(PreyNumToPlotVec(kPrey)).JammingnTimes);
         plot(h1,PREY(PreyNumToPlotVec(kPrey)).xPrey(ixJam), PREY(PreyNumToPlotVec(kPrey)).yPrey(ixJam), ...
             'Color', colStruct.masking , 'Marker', 'x', ...
             'MarkerSize',9, 'LineStyle', 'none' , 'DisplayName', 'PreyClicks' ) % the position of the Prey

     end % for kPrey
 end % if DataToPlot.JamPreyPosFlag

% % % % The Catches of the Prey positions
% % %     if DataToPlot.CatchesPosFlag
% % %         CatchPreyTimes = DataToAnalyze.BAT(BatToPlot).CatchPreyTimes./SampleTime;
% % %         CatchPreyPos = DataToAnalyze.BAT(BatToPlot).CatchPreyPos;
% % %         
% % %         if ~isempty(CatchPreyPos)
% % %             for nCatch = 1: DataToAnalyze.BAT(BatToPlot).NumberOfCatches
% % %                 Index = find(IDn == round(CatchPreyTimes(nCatch)),1); % Check if the catch in ZOOM
% % %                 if ~isempty(Index) 
% % %                     plot(h1, CatchPreyPos(nCatch,1), CatchPreyPos(nCatch,2), 'k+','MarkerSize',16,'linewidth',4)
% % %                 end % if ~isempty(Index)
% % %             end % for nCatch
% % %         end % if ~isempty(CatchPreyPos)
% % %     end % DataToPlot.ManPosFlag

% axis(h1,[ -1,TerrainParams.Xmax+1,...
%         -1, TerrainParams.Ymax+1]);

all_xbats = [DataToAnalyze.BAT.xBatPos];
all_ybats = [DataToAnalyze.BAT.yBatPos];
Xmin = min([-1, all_xbats]); Xmax = max([TerrainParams.Xmax, all_xbats]); 
Ymin = min([-1, all_ybats]); Ymax = max([TerrainParams.Ymax, all_ybats]); 

axis(h1,[ Xmin-0.5,Xmax+0.5,...
        Ymin-0.5, Ymax+0.5]);
if strcmp( TerrainParams.Type, 'Room-L obs')
    xlim([-0.2 10.2]);
    ylim([-0.2 9.2]);
    title(['Cave Exit Example', num2str(TotalBatsNumber), 'Bats'])
end

    