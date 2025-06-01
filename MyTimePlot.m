
function [] = MyTimePlot(HTimePlot,DataToAnalyze,DataToPlot,...
    SimParams, TerrainParams, Terrain, ZoomFlag, ZoomIDn)

h2 = HTimePlot;
hold(h2,"on")
SampleTime = DataToAnalyze.AllParams.SimParams.SampleTime;
xyResolution = DataToAnalyze.AllParams.SimParams.xyResolution;
MaxTime = DataToAnalyze.AllParams.SimParams.SimulationTime;
TimeVec = [0:SampleTime:MaxTime]; % Time Vector
ColorByManFlag = 1;

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

StartTimeIndex = 1;
EndTimeIndex = MaxTime/SampleTime+1;
%%% Dealing with Zoom %%%
if DataToPlot.TimeFilterFlag && ~ZoomFlag % 
    cla(h2);
    StartTimeIndex = max(DataToPlot.StartTimeFilter./SampleTime , 1);
    EndTimeIndex = DataToPlot.EndTimeFilter./SampleTime+1;
end % if DataToPlot.TimeFilterFlag && ~ZoomFlag
IDn = StartTimeIndex:EndTimeIndex;
if ZoomFlag
    cla(h2);
    IDn = ZoomIDn;
end % ZoomFlag
TimeVec = TimeVec(IDn);

%%% Setting Colors
TotalBatsNumber = DataToAnalyze.AllParams.SimParams.TotalBatsNumber;
TotalPreysNumber = DataToAnalyze.AllParams.SimParams.TotalPreysNumber;

[colStruct] = myColorScheme(SimParams.TestMode, TotalBatsNumber, TotalPreysNumber);
ColorOfPreys     = colStruct.PreysPos;
ColorHuntingMan  = colStruct.Approach; % [ 0 1 0];  % Green
ColorObsMan      = colStruct.ObsMan; % [0 0.7 0]; % Green
ColorFindingPrey = colStruct.detections; %[1 0 1]; % magenta
ColorForagingMan = colStruct.ForagingMan; % [0 0 0]; %Black


%%% The Plots %%%
for kBat = 1:NumberOfBatsToPlot
    BatToPlot = BatNumToPlotVec(kBat);
    if  DataToPlot.FlightVelocity
        %     Direction = atan(diff([xBatPos,xBatPos(end)])./diff([yBatPos,yBatPos(end)]));
        BatVelocity = DataToAnalyze.BAT(BatToPlot).BatVelocity * xyResolution / SampleTime; % m/sec
        hold(h2,"on")
        plot(h2,TimeVec,BatVelocity, 'b-'  )
        axis auto
        %     axis(h2,[TimeVec(1) ,TimeVec(end),min(DataToAnalyze.BatDirection)-0.1,...
%         max(DataToAnalyze.BatDirection)+0.1]);
    end %  DataToPlot.FlightDirection
    
    if  DataToPlot.SonarTimes
        if ColorByManFlag
            MyTimeColorByAnyFieldPlot( DataToAnalyze.BAT(BatToPlot).SonarPulses(IDn), 'FlightStageCell',DataToAnalyze.BAT(BatToPlot) , DataToAnalyze.AllParams, h2)
%             MyTimeColorByManueverPlot( DataToAnalyze.BAT(BatToPlot).SonarPulses(IDn), DataToAnalyze.BAT(BatToPlot) , DataToAnalyze.AllParams, h2)
        else %if ColorByManFlag
             plot(h2,TimeVec ,  DataToAnalyze.BAT(BatToPlot).SonarPulses(IDn), colStruct.EcholocatePos)
        end %  if ColorByManFlag
        
        axis('auto');
    end %  DataToPlot.SonarTimes
    
   
    
%     if  DataToPlot.RadialAccel
%          plot(h2,TimeVec ,  DataToAnalyze.BAT(BatToPlot).BatDteta,'r' )
%     end
    
%     if  DataToPlot.ObsManuevers
%         bar(h2,TimeVec ,  DataToAnalyze.BAT(BatToPlot).ObsManueverFlag(IDn), 'm' )
%         axis(h2,[ TimeVec(1) ,TimeVec(end),-0.1, 1.1]);
%     end
    % PRFVec
    if DataToPlot.PRFVec
%         plot(h2,TimeVec ,  DataToAnalyze.BAT(BatToPlot).PRFVec(IDn), 'g' )
        
        IPIVec= zeros(size(DataToAnalyze.BAT(BatToPlot).PRFVec));
        [~, PulseInd, Prfs] = find(DataToAnalyze.BAT(BatToPlot).PRFVec);
        IPIVec(PulseInd) = 1e3./Prfs;
        if ColorByManFlag
%             MyTimeColorByManueverPlot( DataToAnalyze.BAT(BatToPlot).PRFVec, DataToAnalyze.BAT(BatToPlot) , DataToAnalyze.AllParams, h2)
%             MyTimeColorByManueverPlot( IPIVec, DataToAnalyze.BAT(BatToPlot) , DataToAnalyze.AllParams, h2)
            MyTimeColorByAnyFieldPlot( IPIVec, 'FlightStageCell', DataToAnalyze.BAT(BatToPlot) , DataToAnalyze.AllParams, h2)
        else %if ColorByManFlag
             plot(h2,TimeVec ,  DataToAnalyze.BAT(BatToPlot).PRFVec(IDn), colStruct.ownCalls )
        end %  if ColorByManFlag
        
    end % if  DataToPlot.PRFVec
    
    %  Plot Pulse Duration
    if DataToPlot.PulseDuration
        PulseDurationMS =  DataToAnalyze.BAT(BatToPlot).PulseWidthVec.*SampleTime*1000;
%         MyTimeColorByManueverPlot( PulseDurationMS, DataToAnalyze.BAT(BatToPlot) , DataToAnalyze.AllParams, h2)
        MyTimeColorByAnyFieldPlot( PulseDurationMS, 'FlightStageCell', DataToAnalyze.BAT(BatToPlot) , DataToAnalyze.AllParams, h2)
    end % if DataToPlot.PulseDuration
    
    
    % Distance and angle from Prey
    if  DataToPlot.Bat2PreyDistFlag
        for kPrey = 1:NumberOfPreysToPlot
            PreyToPlot = PreyNumToPlotVec(kPrey);
            % Check wether the prey was caught and cut the e-relevant data
            IDn1 = IDn;
%             if DataToAnalyze.PREY(kPrey).IsCaught
%                 nCaughtTime = DataToAnalyze.PREY(kPrey).CaughtTime / SampleTime ;
%                 
%                 if (nCaughtTime) < IDn(1)
%                     IDn1 =[];
%                 else % if CaughtTime < IDn(1)
%                        IDn1 = uint32(IDn(1): nCaughtTime);   
%                 end % CaughtTime < IDn(1)
%                 
%             end %if DataToAnalyze.PREY(kPrey).IsCaught
            
            TimeVecXX = TimeVec(IDn1);  
            Dist2Prey = DataToAnalyze.BAT(BatToPlot).Vector2Preys(PreyToPlot).Dist2Prey(IDn1);
            plot(h2,TimeVecXX , Dist2Prey , 'color', ColorOfPreys(PreyToPlot,:), 'linewidth', 3 )
            hold on
%             if ColorByManFlag
% %                 MyTimeColorByManueverPlot( Dist2Prey, DataToAnalyze.BAT(BatToPlot) , DataToAnalyze.AllParams, h2,...
% %                 'Marker', '.', 'LineStyle','none','MarkerSize',4)
%                  MyTimeColorByAnyFieldPlot( Dist2Prey, 'FlightStageCell' ,DataToAnalyze.BAT(BatToPlot) , DataToAnalyze.AllParams, h2,...
%                 'Marker', '.', 'LineStyle','none','MarkerSize',4)
%             else %if ColorByManFlag
% %                plot(h2,TimeVecXX , Dist2Prey) 
%             end %  if ColorByManFlag
        end % for
        axis(h2, 'auto')
    end % if
    
    
    
    % Angle to Prey
    if  DataToPlot.Bat2PreyAngleFlag
        for kPrey = 1:NumberOfPreysToPlot
            PreyToPlot = PreyNumToPlotVec(kPrey);
            Angle2Prey = DataToAnalyze.BAT(BatToPlot).Vector2Preys(PreyToPlot).Angle2Prey(IDn);
            plot(h2,TimeVec , Angle2Prey , 'k' )
        end % for
        axis auto
    end
    
   
    
    % The Velocity of the Bat
    if DataToPlot.FlightVelocity
       plot(h2,TimeVec ,DataToAnalyze.BAT(BatToPlot).BatVelocity);
    end
    
    
    
    
     % AllPulses plots
%     L = min(length(TimeVec), length(DataToAnalyze.PreyEchosVec));
    plotDetectionFlg = 1;

    if DataToPlot.BatAllPulsesFlag
        plot(h2, TimeVec, 10*log10(abs(DataToAnalyze.BAT(BatToPlot).BatSonarEchosMat(1,IDn) )), 'color', colStruct.ownCalls, 'linewidth',1)
        axis auto
    end
    
    % The Interfernce
%     DataToPlot.InterferencePlotflag = 0;
    if DataToPlot.InterferencePlotFlag
        DetectionTH = DataToAnalyze.AllParams.BatSonarParams.PulseDetectionTH;
        DetectionTHVec = DetectionTH*ones(1,length(TimeVec));
        
        plot(h2,TimeVec , 10*log10(abs(DataToAnalyze.BAT(BatToPlot).AllInterPulses)), '-', 'color', colStruct.masking, 'linewidth',1 )
        plot(h2, TimeVec, DetectionTHVec, 'g-','linewidth',2);
        
    end %if DataToPlot.InterferencePlotflag
    
    % Obs Echos plots
    if DataToPlot.ObsEchosFlag
        plot(h2,TimeVec , 10*log10(abs(DataToAnalyze.BAT(BatToPlot).ObsEchosVec(IDn))), 'o-', 'color', colStruct.clutter, 'MarkerSize',2 )
        axis auto
        if plotDetectionFlg && DataToAnalyze.AllParams.SimParams.DetectObs
            detTimes = [DataToAnalyze.BAT(BatToPlot).ObsFindsStruct.DetectedTimes];
            ix       = round(detTimes);
            plot(h2, detTimes*SampleTime, 10*log10(DataToAnalyze.BAT(BatToPlot).ObsEchosVec(ix)),'*', 'color', colStruct.detections, 'LineWidth', 2)
        end
    end
    
    % Prey Echos plots
    if DataToPlot.PreyEchosFlag
        plot(h2,TimeVec , 10*log10(abs(DataToAnalyze.BAT(BatToPlot).PreyEchosVec(IDn))), '+', 'color', colStruct.Prey,'MarkerSize',4 )
        axis auto
        if plotDetectionFlg && DataToAnalyze.AllParams.SimParams.DetectPrey
            detTimes = [DataToAnalyze.BAT(BatToPlot).PreyFindsStruct.DetectedTimes];
            ix       = round(detTimes);
            plot(h2, detTimes*SampleTime, 10*log10(DataToAnalyze.BAT(BatToPlot).PreyEchosVec(ix)),'*', ...
                'color', colStruct.detections, 'LineWidth', 2)
        end
    end

     % Consps Echos plots
    if DataToPlot.ConspsEchosFlag
        plot(h2,TimeVec , 10*log10(abs(DataToAnalyze.BAT(BatToPlot).Consps_EchosVec(IDn))), '.-', 'color', colStruct.conspsEchoes,'MarkerSize',3 )
        axis auto
        if plotDetectionFlg && DataToAnalyze.AllParams.SimParams.DetectConsps
            detTimes = [DataToAnalyze.BAT(BatToPlot).Consps_FindsStruct.DetectedTimes];
            ix       = round(detTimes);
            plot(h2, detTimes*SampleTime, 10*log10(DataToAnalyze.BAT(BatToPlot).Consps_EchosVec(ix)),'*', 'color', colStruct.detections, 'LineWidth', 2)
        end
    end
    

    %%% Add Transmittigng Pulse num 
    if DataToPlot.PreyEchosFlag || DataToPlot.ObsEchosFlag || DataToPlot.InterferencePlotFlag || DataToPlot.BatAllPulsesFlag || DataToPlot.FlightVelocity
        TxTimes = [DataToAnalyze.BAT(BatToPlot).TransmittedPulsesStruct.StartPulseTime] * SampleTime;
%         t0 = text(-0.085*max(TimeVec), -5, 'time:','FontSize', 8);
        t1 = text(TxTimes, -10*ones(size(TxTimes)), string([DataToAnalyze.BAT(BatToPlot).TransmittedPulsesStruct.PulseNum]),'FontSize', 8);
        t2 = text(-0.15*max(TimeVec), -10, 'PulseNum:','FontSize', 8);

    end % if DataToPlot.PreyEchosFlag || DataToPlot.ObsEchosFlag || DataToPlot.InterferencePlotFlag || BatAllPulsesFlag

     % Bat Manuvers Plot
    if  DataToPlot.HuntingManuevers
        
        % the Indices Of the Manuevers
        RefPlot= zeros(1,length(TimeVec));
        MyTimeColorByAnyFieldPlot( RefPlot, 'FlightStageCell', DataToAnalyze.BAT(BatToPlot), DataToAnalyze.AllParams, h2,...
            'Marker', 'o', 'LineStyle','none','MarkerSize',6 );
%         MyTimeColorByManueverPlot( RefPlot, DataToAnalyze.BAT(BatToPlot) , DataToAnalyze.AllParams, h2,...
%             'Marker', 'o', 'LineStyle','none','MarkerSize',6 );
        
    end %  DataToPlot.HuntingManuevers
    
    
    % The Times Of Jamming
    DataToPlot.JammingTimesFlag = 1;
    if DataToPlot.JammingTimesFlag 
        switch DataToPlot.DetectObjToAnalyze
            case 'Prey'
                JammedTimes = double(DataToAnalyze.BAT(BatToPlot).InterReportStrctOnLine.TotalInterferenceTimes); % ; .*SampleTime ;
            case 'Conspecific'
                JammedTimes = double(DataToAnalyze.BAT(BatToPlot).InterReportStrctOnLine.conspsTotalInterferenceTimes); % ; .*SampleTime ;
            case 'Obs'
                 JammedTimes = double(DataToAnalyze.BAT(BatToPlot).InterReportStrctOnLine.obsTotalInterferenceTimes   ); % ; .*SampleTime ;
        end
        if ~isempty(JammedTimes)
            plot(h2, JammedTimes, 40*ones(size(JammedTimes)),'db','MarkerSize',8,'color', colStruct.masking,'linewidth',1.5)
        end % if ~isempty(JammedTimes)
    end % if DataToPlot.JammingTimesFlag
    
    % Catching Times
    DataToPlot.CatchTimesFlag = 1;
    if DataToPlot.CatchTimesFlag
        CatchPreyTimes = DataToAnalyze.BAT(BatToPlot).CatchPreyTimes;
        if ~isempty(CatchPreyTimes)
            plot(h2, CatchPreyTimes, zeros(size(CatchPreyTimes)),'+', 'color', colStruct.catchPrey, 'MarkerSize',10,'linewidth',4);
        end % if ~isempty(JammedTimes)
    end % if DataToPlot.CatchTimesFlag
    
end % for BatNum = 1:NumberOfBatsToPlot

