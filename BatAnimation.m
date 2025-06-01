function[]  = BatAnimation(varargin)
%BatAnimation(DataToAnalyze, , SimParams, DataToPlot AnimateSpeed)
%This function will animate the bat and prey movemet in time
%
%Inputs:
%   DataToAnalyze: A struct with the vectors of the bat and the prey
%   movements:
%      BatPos   BatDirection    BatDirectionChange  SonarTimesVec
%      ObsManueverFlag      HuntingFlag     ObsticlesFound
%      PreyFound    PRFVec  PulseWidthVec   PreyPos     PreyDirection
%      CatchPreyPos     Dist2Prey       Angle2Prey
%   DataToPlot: A struct with flag which data to plot
%      PreyPositionFlag TerrainFlag BatPosFlag  SonarPosFlag    ManPosFlag
%      TerrainFindFlag  EndTimeFilter   StartTimeFilter     TimeFilterFlag
%      CatchesPosFlag   HuntingManuevers    FlightDirection SonarTimes
%      ObsManuevers     SonarPulses     PRFVec      Bat2PreyDistFlag
%      Bat2PreyAngleFlag    PreyFindsFlag
%   SimParams: A struct with genearal parameters of the simulation 
%   AnimateSpeed: the time speed, the higher the faster 
%
%   Created 10/10/17 by Omer

NumOfInputs = nargin;
%
% default Parameters
SampleTime = 1e-4;
xyResolution = 0.01;
MaxTime = 10/SampleTime;

AnimateSpeed = 5000;
TraceAnimationSize = 5*AnimateSpeed;
hFigure = figure();
%subplot(2,1,1)
axis( [-5,15,-5, 15]);
title('Bat Animation');
% % subplot(2,1,2)
% % title('time animation')
% % axis( [-0.1, 10, 0, MaxTime]);
if NumOfInputs >= 1
        DataToAnalyze = varargin{1};
        xBatPos = DataToAnalyze.BAT(1).xBatPos;
        yBatPos = DataToAnalyze.BAT(1).yBatPos;
%         xPreyPos = DataToAnalyze.PreyPos(1,:);
%         yPreyPos = DataToAnalyze.PreyPos(2,:);
        SampleTime =  DataToAnalyze.AllParams.SimParams.SampleTime;
        MaxTime = DataToAnalyze.AllParams.SimParams.SimulationTime./SampleTime;
        NumberOfPreys= DataToAnalyze.AllParams.SimParams.TotalPreysNumber;
        xyResolution = DataToAnalyze.AllParams.SimParams.xyResolution;
        PreyPositionFlag = 1;
        BatPositionFlag = 1;
        ManueversFlag = 0;
%         Dist2PreyFlag = 1;
% % %         if NumOfInputs == 2
% % %             % SimParams
% % %             SimParams = varargin{2};
% % %             NumberOfPreys= SimParams.TotalPreysNumber;
% % %             SampleTime =  SimParams.SampleTime;
% % % %             MaxTime =  SimParams.SimulationTime/SampleTime;
% % %             xyResolution = SimParams.xyResolution;
% % %         end
end % if NofInputs >= 1
        
  
    
hBatPos = animatedline('Marker','o', 'color','k');

Color = ['b','c','g', 'r', 'k', 'y', 'm'];
for nPrey = 1:NumberOfPreys
   
    hPreyPos{nPrey} = animatedline('Marker','+', 'color',Color(nPrey));
end
if ManueversFlag
    hHuntingPos = animatedline('Marker','o', 'color','r');
    hObsManPos =  animatedline('Marker','.', 'color','r');
    HuntingFlag = DataToAnalyze.BAT(1).HuntingFlag;
    ObsManueverFlag = DataToAnalyze.BAT(1).ObsManueverFlag;
end % if HuntingManueversFlag

% % % hDist2PreyTime = animatedline('Marker','.', 'color','r');
% % % Dist2Prey = DataToAnalyze.BAT(1).Vector2Preys.Dist2Prey;

PREY(NumberOfPreys) = struct('xPreyVec',[], 'yPreyVec',[]);
for tt = 1: AnimateSpeed: MaxTime-AnimateSpeed
    timeToPlot= max(1, tt-TraceAnimationSize);
    %
% %     subplot(2,1,1) % X_Y animation
    if BatPositionFlag
        xBatVec = xBatPos(tt:tt+AnimateSpeed-1);
        yBatVec = yBatPos(tt:tt+AnimateSpeed-1);
        addpoints(hBatPos, xBatVec, yBatVec)
    end % if BatPositionFlag
    if PreyPositionFlag
        for nPrey = 1:NumberOfPreys
            PREY(nPrey).xPreyVec = DataToAnalyze.PREY(nPrey).xPrey(timeToPlot:tt+AnimateSpeed-1);
            PREY(nPrey).yPreyVec = DataToAnalyze.PREY(nPrey).yPrey(timeToPlot:tt+AnimateSpeed-1);
            addpoints(hPreyPos{nPrey}, PREY(nPrey).xPreyVec, PREY(nPrey).yPreyVec);
        end
    end % if BatPositionFlag
    if ManueversFlag
        HuntingFlagVec = HuntingFlag(timeToPlot:tt+AnimateSpeed-1);
        ObsManueverVec = ObsManueverFlag(timeToPlot:tt+AnimateSpeed-1);
        xHunting = xBatPos(timeToPlot + find(HuntingFlagVec));
        yHunting = yBatPos(timeToPlot + find(HuntingFlagVec));
        xObsManuever = xBatPos(timeToPlot + find(ObsManueverVec));
        yObsManuever = yBatPos(timeToPlot + find(ObsManueverVec));
        addpoints(hHuntingPos, xHunting, yHunting)
        addpoints(hObsManPos, xObsManuever, yObsManuever);
    end % if HuntingManueversFlag
    drawnow limitrate
    clearpoints(hBatPos)
%     clearpoints(hPreyPos)
    %
%     subplot(2,1,2)  % Time Vec Animation
%     if Dist2Prey
%         TimeVec = [tt:tt+AnimateSpeed-1];
%         Dist2PreyVec = Dist2Prey(TimeVec);
%         addpoints(hDist2PreyTime, Dist2PreyVec, TimeVec)
%     end % if Dist2Prey
%     drawnow %limitrate
    %
end % for tt 

% Final Plot
hold ('on')
plot(xBatPos, yBatPos,'ok')
% plot(xPreyPos, yPreyPos,'+c')
% % % plot(xBatPos(find(HuntingFlag)),yBatPos(find(HuntingFlag)),'or')
% plot(xBatPos(find(ObsManueverFlag)),yBatPos(find(ObsManueverFlag)),'or')