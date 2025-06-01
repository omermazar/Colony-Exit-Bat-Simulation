function [] = MyBatUmweltOnLine( BAT , AllParams, varargin)

% this function will build and plot the umwelt of the bats
% Inputs:
%   BAT = BatDATA.BAT(BatNum)
%   AllParams - struct with the required data fro the simularion
%   optional {varagin}- DataToPlot = 1 - 'Distance' 2- 'RxPower' 3-'SIR' 4- 'Angle2Prey' 
%                                   5- 'DiffAngle'
%    
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
%%
NofInputs = nargin; 
% 
if NofInputs ==3
    DataToPlotCode = varargin{1};
else % if NofInputs ==3
    DataToPlotCode= 1;
end % if NofInputs ==3

    
InterReportStrct = BAT.InterReportStrctOnLine;
EchosFromObsStruct  = BAT.EchosFromObsStruct;

PreyFindsStruct = BAT.PreyFindsStruct;
if BAT.NumberOfCatches > 0
    CatchPreyNum = BAT.CatchPreyNum;
    CatchPulseNum = BAT.CatchPulseNum;
end % if BAT.NumberOfCatches


%%


%%% The Obtacles Detections %%%
%%
TotalNumberOfPulses = InterReportStrct.TotalNumberOfPulses;
ObsDetecitionsVec = zeros(1,TotalNumberOfPulses);
ObsValue = AllParams.SimParams.TotalPreysNumber + 0.5;  % the obs found should be drawn after the preys
if ~isempty (EchosFromObsStruct)
    PulseOfDetecctedObs = [EchosFromObsStruct.TransmittedPulseNum];
    ObsDetecitionsVec(PulseOfDetecctedObs) = ObsValue ;
end  % if ~isempty (EchosFromObsStruct)

%%

%%% Matrix of Detects OnLine %%%
%%
%%% Init variables 
SampleTime = AllParams.SimParams.SampleTime;
xyResolution = AllParams.SimParams.xyResolution;

Dist2PreyfullVec = [PreyFindsStruct.Dist2DetectedPrey];
HuntedPreyVec = zeros(1,TotalNumberOfPulses);
NumOfPreys = AllParams.SimParams.TotalPreysNumber;
MatrixOfDetects = zeros(NumOfPreys, TotalNumberOfPulses);
MatrixOfMasking = zeros(NumOfPreys, TotalNumberOfPulses);
SIRMatrix = zeros(NumOfPreys, TotalNumberOfPulses);
Dist2PreyMatrix = -(max(Dist2PreyfullVec)+100)*ones(NumOfPreys, TotalNumberOfPulses);
RxPowerMatrix = zeros(NumOfPreys, TotalNumberOfPulses);
Angle2PreyMatrix = -180*ones(NumOfPreys, TotalNumberOfPulses);
AngularVel2PreyMatrix = zeros(NumOfPreys, TotalNumberOfPulses);
RelativeVelocity2PreyMatrix = zeros(NumOfPreys, TotalNumberOfPulses);
IPIsecVector = zeros(1, TotalNumberOfPulses);

%%% Matrices of all Distances and angles at the start of each pulse
RelPositionStruct = BAT.PreysRelativePositionVecStr;  

AllDistances = [RelPositionStruct.Distances];
AllDistancesMat = reshape(AllDistances, NumOfPreys, TotalNumberOfPulses);
% adjust to the end of pulse
AllDistancesMat = [ AllDistancesMat(:,2:end), 5*ones(NumOfPreys,1)];
DiffAllDistances = diff(AllDistancesMat')';
DiffAllDistances(:,end+1) = zeros(NumOfPreys,1);

AllAngles = [RelPositionStruct.Bat2TargetRelativeAngle];
AllAnglesMat = reshape(AllAngles, NumOfPreys, TotalNumberOfPulses);
AllAnglesMat = [ AllAnglesMat(:,2:end), 5*ones(NumOfPreys,1)];
DiffAllAngles = diff(AllAnglesMat')';
DiffAllAngles(:,end+1) = zeros(NumOfPreys,1);

% fullfill the results
for k = 1:TotalNumberOfPulses
    IPIsecVector(k) = BAT.TransmittedPulsesStruct(k).IPItoNextPulse * SampleTime;
    DetectedPreys = PreyFindsStruct(k).DetecectedPreyWithOutInterference; 
    MissedPreysByProbabilty = PreyFindsStruct(k).MissedPreysByProbabilty;
    
    if ~isempty(DetectedPreys)
        if~isempty(PreyFindsStruct(k).PreyNumToHunt)
            HuntedPreyVec(k) = PreyFindsStruct(k).PreyNumToHunt;
        end % if~isempty(PreyFindsStruct(k).PreyNumToHunt)
        MatrixOfDetects(DetectedPreys,k)= 1;
        MatrixOfMasking(PreyFindsStruct(k).MaskedPreys,k) = 1;
        SIRMatrix(DetectedPreys,k) = PreyFindsStruct(k).SIROfDetectedPreys;
        SIRMatrix( SIRMatrix == 0 ) = -60;
        
        Dist2PreyMatrix(DetectedPreys,k) = -PreyFindsStruct(k).Dist2DetectedPrey;
        
        Angle2PreyMatrix(DetectedPreys,k) = -abs(PreyFindsStruct(k).Angle2DetectedPrey)*180/pi;
        RxPowerMatrix(DetectedPreys,k) = PreyFindsStruct(k).RxPowerOfDetectedPreys;
        
        
        % the relative veloctuy and angular accelaration
        RelativeVelocity2PreyMatrix(DetectedPreys,k) = ...
            DiffAllDistances(DetectedPreys,k) * xyResolution ./ IPIsecVector(k)  ; 
        
        AngularVel2PreyMatrix(DetectedPreys,k) = ...
            abs(DiffAllAngles(DetectedPreys,k)) *(180/pi) ./ IPIsecVector(k) ;
        
    end % if ~isempty(DetectedPreys)
    if ~isempty(MissedPreysByProbabilty)
        MatrixOfDetects(MissedPreysByProbabilty,k)= 0.3;
    end % if ~isempty(DetectedPreys)
end %for k
 
% DiffAngle2PreyMatrix = diff(Angle2PreyMatrix')' * 180/pi ./
% IPImsecVector; Error

MaskPlotConst = 0.7;
DetectionJamMat = MatrixOfDetects- MaskPlotConst*MatrixOfMasking;
[MaskedPreyVec, MaskedPulseVec ] = find(MatrixOfMasking);


%%
%%% Stage of Manuver
FlagStageToPlot = 0;

if FlagStageToPlot
    uStages = BAT.InterReportStrctOnLine.SumStagesStruct.Stage;
    for kk= 1:length(uStages)
        StageInd.(uStages{kk}) = ...
            cell2mat(BAT.InterReportStrctOnLine.SumStagesStruct.StagePulsesCell(kk, : ));
        PreyByStage.(uStages{kk}) = HuntedPreyVec(StageInd.(uStages{kk}));  
    end % for kk
end % if FlagStageToPlot



%%

%%% PLOTS 
%%


hF = figure; % umwelt by function
hF.Position = [520 100 1000 500];

hold on
% 'Distance' 'RxPower' 'SIR' 'Angle2Prey' 
switch DataToPlotCode
    case 1 % 'Distance'
        Th = title('distance from prey [cm]');
        
        imagesc(Dist2PreyMatrix); colorbar
        caxis([-400 , 0]);
%         caxis([min(Dist2PreyMatrix(:)) , 0]);
    
    case 2 % 'RxPower'
        Th = title('ReceivedPower from prey [dB]');
        imagesc(RxPowerMatrix); colorbar
        caxis([0, max(RxPowerMatrix(:))]);
   
    case 3 %'SIR'
        Th = title('Signal to Interference [dB]');
        imagesc(SIRMatrix); colorbar
        caxis([min(SIRMatrix(:)) , max(SIRMatrix(:))]);
    
    case 4 %'Angle2Prey'
        Th = title('Angle to prey [deg]');
        imagesc(Angle2PreyMatrix); colorbar
        caxis([ min(Angle2PreyMatrix(:)), 0]);
    
    case 5 %'DiffAngle'
        Th = title('Angular Velocity to prey (Manuever Intensity) [deg/sec]');
        imagesc(AngularVel2PreyMatrix); colorbar
        caxis([ 0, 150]);
    
    case 6 %'Relative Velocity'
        Th = title('Relative Velocity to prey [m/sec]');
        imagesc(RelativeVelocity2PreyMatrix); colorbar
        caxis([ -3, 3]);    
end % switch DataToPlot

% moving up the title
% Th_pos= get(Th, 'position');
% set(Th, 'position', get(Th, 'position') + [ 0 0.2 0]);

hold on
ylabel('prey number') ; xlabel('Pulse Number')
plot(MaskedPulseVec, MaskedPreyVec, 'd', 'MarkerSize',8,'LineWidth',6,'color','r')

plot(HuntedPreyVec,'*g','MarkerSize',8,'LineWidth',6)
plot(ObsDetecitionsVec,'*w', 'MarkerSize',8,'LineWidth',6)

Ltext1 = {'Masked Echo', 'Hunted Prey' ,'Obstacle Detections '};
L = legend(Ltext1); 
% L = legend('Masked Echo', 'Hunted Prey' ,'Obstacle Detections '); 

% Catches
if BAT.NumberOfCatches
    plot( CatchPulseNum, CatchPreyNum, '+', 'MarkerSize',10,'LineWidth',6,'color','k')
    Ltext2 = {Ltext1{:}, 'Catches'};
    L = legend(Ltext2); 
%     L = legend('Masked Echo', 'Hunted Prey' ,'Obstacle Detections ' ,'Catches'); 
%     CatchTxtLeg= 'Catches';
else % if BAT.NumberOfCatches
    L.String = {L.String{1:end-1}};
%     CatchTxtLeg=  '';
end % if BAT.NumberOfCatches

%%% Stages Plot
if FlagStageToPlot
    ColorOfStages = {'g', 'm', 'b' , 'c', 'y'};
    for kk= 1:length(uStages)
        plot(ax1, StageInd.(uStages{kk}), PreyByStage.(uStages{kk}),'*','MarkerSize',8, 'color', ColorOfStages{kk}) 
    end % for kk
    L = legend(Ltext2{:},(uStages{:}));
end % if FlagStageToPlot

axis( [ 0, TotalNumberOfPulses, 0, ObsValue+0.5 ] );
% L = legend('Masked Echo', 'Hunted Prey' ,'Obstacle Detections ' ,'Catches'); 
L.Location = 'northeastoutside';
text(0,ObsValue,'Obstacle', 'HorizontalAlignment','right')
%%

%%% Add Time SCale on X-axis up
% x1 = [BAT.TransmittedPulsesStruct.PulseNum]; % the lower x-axis by pulse num
% Pt = [BAT.TransmittedPulsesStruct.StartPulseTime]; % the upper shoul be by pulse time
% 
% xt1 = get(gca, 'XTick'); % the positions of the x-axis
% Ind = ismember(x1, xt1);
% Ptx1 = Pt(Ind);
% Ptx1= [0, Ptx1, Pt(end)] * SampleTime;
% 
% ax1 = gca;
% ax1_pos = ax1.Position;% position of first axes
% 
% % sets the upper x-axis
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% xt2 =get(ax2, 'XTick');
% xtup = linspace(xt2(1), xt2(end), length(Ptx1));
% 
% set(gca, 'XTick', xtup, 'XTickLabel', Ptx1)

%%

%%% PushButtons
%%
hPopMenu = uicontrol('Style', 'popupmenu' , ...
    'String', {'Distance'; 'RxPower'; 'SIR'; 'Angle2Prey'; 'AngularVelocity' ; 'RelativeVelocity'}, ...
    'Position', [780 10 80 20], ...
    'Callback', @popup_menu_Callback);
% PopMenu.Callback = 'MyBatUmweltOnLine( BatDATA.BAT(1) , AllParams, 3 )';


%%% Pop-up menu callback. Read the pop-up menu Value property to
%  determine data to plot
    function popup_menu_Callback(source,eventdata)
        % Determine the selected data set.
        Str = get(source, 'String');
        PopupValue = get(source,'Value');
%         popupValue= PopupContents{get(source,'Value')};
        MyBatUmweltOnLine( BAT , AllParams, PopupValue )
    end % function popup_menu_Callback(source,eventdata)
%%
% plot Stages Button
hPlotStages = uicontrol('Style', 'pushbutton' ,'String', 'plot stages', ...
    'Position', [780 50 80 30], ...
    'Callback', @PlotStages_Callback);


%%% PushButton callback. Read the pop-up menu Value property to
%  determine data to plot
    function PlotStages_Callback(source,eventdata)
        % Determine the selected data set.
        ColorOfStages = {'k', 'r', 'b' ,[0.5 0.5 0],  [0.5 0.0 0.5]};
        textYpos  = [4:-0.5:-1];
        uStages = BAT.InterReportStrctOnLine.SumStagesStruct.Stage;
        for kk= 1:length(uStages)
            StageInd.(uStages{kk}) = ...
                cell2mat(BAT.InterReportStrctOnLine.SumStagesStruct.StagePulsesCell(kk, : ));
            PreyByStage.(uStages{kk}) = HuntedPreyVec(StageInd.(uStages{kk}));
            
            plot(StageInd.(uStages{kk}), PreyByStage.(uStages{kk}),...
                '*','MarkerSize',8, 'color', ColorOfStages{kk})
            t = text(-5,textYpos(kk), uStages{kk}, 'color', ColorOfStages{kk},...
                'HorizontalAlignment','right');
        end % for kk
%         L = legend((uStages{:}));
%         L.Location = 'northeastoutside';
    end % function popup_menu_Callback(source,eventdata)

%%
end % main function
% Online detection plot
% figure
% hold on
% imagesc(MatrixDetect)
% plot(SearchPulsesNum,PreyToHunt(SearchPulsesNum),'o','MarkerSize',8,'LineWidth',4,'color',[1 1 0])
% plot(AproachPulsesNum,PreyToHunt(AproachPulsesNum),'o','MarkerSize',8,'LineWidth',4,'color',[ 0 1 0])
% end

