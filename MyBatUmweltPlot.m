function [DetectionJamMat] = MyBatUmweltPlot( BAT )

% this function will build and plot the umwelt of the bats
% outputs - 
%    1) matrix of detections and Interdernce of preys and Obstacles - DetectionJamMat
%           X - Pulse num 
%           Y(1IN) - Prey Num, Y(end)- Obstacle detection
%           Mat(Npulse, Prey) = Detection - 0.7*Interference\
%           ObstacleValue= detection
%    2) matrix of stages - ObsMan or Appraoch
%             Mat(Npulse, Prey)  = 1 if Buzz to that Prey 
%                                = 0.5 if Appraoching that Prey
%             Mat(Npulse, Obs)   = 0.5 if ObsMan
%    3) Vector Of Catches = [NPulse, Prey)

%%% Input DATA from BAT
InterReportStrct = BAT.InterReportStrct;
EchosFromObsStruct  = BAT.EchosFromObsStruct;
if BAT.NumberOfCatches > 0
    CatchPreyNum = BAT.CatchPreyNum;
    CatchPulseNum = BAT.CatchPulseNum;
end % if BAT.NumberOfCatches

PreyFindsStruct = BAT.PreyFindsStruct;

%%% The Obtacles Detctions %%%
TotalNumberOfPulses = InterReportStrct.TotalNumberOfPulses;
ObsDetecitionsVec = zeros(1,TotalNumberOfPulses);
PulseOfDetecctedObs = [EchosFromObsStruct.TransmittedPulseNum];

ObsDetecitionsVec(PulseOfDetecctedObs) = 1;

[NumOfPreys, ~] = size(InterReportStrct.PreysDetectionMatrix);
%%% Matrix of Detects OnLine %%%
MatrixDetect = zeros(NumOfPreys, TotalNumberOfPulses);
for k = 1:TotalNumberOfPulses
    DetectedPreys = PreyFindsStruct(k).DetecectedPreyWithOutInterference; 
    MissedPreysByProbabilty = PreyFindsStruct(k).MissedPreysByProbabilty;
    if ~isempty(DetectedPreys)
        MatrixDetect(DetectedPreys,k)= 1;
    end % if ~isempty(DetectedPreys)
    if ~isempty(MissedPreysByProbabilty)
        MatrixDetect(MissedPreysByProbabilty,k)= 0.3;
    end % if ~isempty(DetectedPreys)
end %for k

%%% Matrix Of Detections and Jams Offline analysis%%%
KforJam = 0.7; % parameter to set the color of the Jammed detecions

DetectionJamMat= zeros(NumOfPreys + 1, TotalNumberOfPulses);
DetectionJamMat(1:NumOfPreys , : ) = InterReportStrct.PreysDetectionMatrix - InterReportStrct.PreysInteferenceMatrix*KforJam;
DetectionJamMat(NumOfPreys+1, : ) = ObsDetecitionsVec;

% Vector Of Flight Stages
FlightStegesVec = InterReportStrct.PulsesStageCell;

SearchPulsesNum = find(strcmp('Search',FlightStegesVec));
AproachPulsesNum = find(strcmp('Approach',FlightStegesVec));
ObsManueverPulsesNum = find(strcmp('ObsManuever',FlightStegesVec));

% Vector of the Hunted Prey
PreyToHunt= InterReportStrct.PulseHuntedPreyVec;
% [~,PreyPulseInd,PreyNum] = find(PreyToHunt);

%%% THE PLOT %%%%
%consts
LineSeparateObs= (NumOfPreys+0.5)*ones(1,TotalNumberOfPulses);
PulseNumVec= 1:TotalNumberOfPulses;

figure
title('Bat Decision Matrix')
xlabel('Pulse Number')
ylabel('Prey Num')
Text1 = ' Obstacle';
text(10, NumOfPreys+1, Text1,'color', 'k'); 
hold on

imagesc(DetectionJamMat);
colormap jet;
plot(PulseNumVec,LineSeparateObs,'r', 'LineWidth', 6)

% obsManuver
if ~isempty(ObsManueverPulsesNum)
    plot(ObsManueverPulsesNum, NumOfPreys+1,'o','MarkerSize',8,'LineWidth',4,'color',[0 0.7 0])
end % if ~isempty(ObsManueverPulsesNum)

% Prey - Search and Approach
plot(SearchPulsesNum,PreyToHunt(SearchPulsesNum),'o','MarkerSize',8,'LineWidth',4,'color',[1 1 0])
plot(AproachPulsesNum,PreyToHunt(AproachPulsesNum),'o','MarkerSize',8,'LineWidth',4,'color',[ 0 1 0])
% plot(PreyPulseInd,PreyNum,'+', 'MarkerSize',8)

% Catches
if BAT.NumberOfCatches
    plot( CatchPulseNum, CatchPreyNum, '+',  'MarkerSize',10,'LineWidth',6,'color','k')
end
axis('tight')

% Online detection plot
% figure
% hold on
% imagesc(MatrixDetect)
% plot(SearchPulsesNum,PreyToHunt(SearchPulsesNum),'o','MarkerSize',8,'LineWidth',4,'color',[1 1 0])
% plot(AproachPulsesNum,PreyToHunt(AproachPulsesNum),'o','MarkerSize',8,'LineWidth',4,'color',[ 0 1 0])
% end