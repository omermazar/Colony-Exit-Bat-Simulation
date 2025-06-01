function [] = MyTimeColorByAnyFieldPlot(VectorToPlot, FieldOfValues, BAT, AllParams, FigHandle, varargin)

% plot  the Desired Vector in colors of the manuever
% inputs- VectorToPlot, ManueverTypeCell, SampleTime, FigHandle
%   VectorToPlot- the vector we plot
%   FieldOfValues - The Name of Field (string) of the values we chose the colors from
%       ('ManueverTypeCell' / 'FlightStageCell' ...)
%   BAT = the dat Strucct from the fight
%   AllParmas- the stauct from the simulation - 
%   FigHandle- the handle of the figure
%   VectorType - 'Data' or 'Referenence'
%   optional - varagin - comman to line style, should be in pairs:
%       example: [varargin(1)m,varargin(2)] - [Marker,'+'] 

%%% Getting the Vector of Values by ?It's Nae
%%%% Not General for now
switch FieldOfValues
    case 'ManueverTypeCell'
        ReferenceCell = BAT.ManueverTypeCell;
        Value1 = 'Hunting';
        Value2 = 'ObsMan';
        Value3 = 'Foraging';
        
    case 'FlightStageCell'
%         ReferenceCell = BAT.InterReportStrct.FlightStageCell;
        ReferenceCell =  BAT.InterReportStrctOnLine.FullFlightStages;
        Value1 = 'Approach';
        Value2 = 'ObstacleManuver';
        Value3 = 'Search';
        Value4 = 'Buzz';
        
       uStages = unique(ReferenceCell);
      % remove empty elments in the cell
      uStages =uStages(~cellfun('isempty',uStages));
        
end % switch FieldOfValues



NumOfExtraInputs = length(varargin);
NumOfComms = floor(NumOfExtraInputs/2);
PlotCommand = repmat({''},1,NumOfComms);
PlotComVal =  repmat({''},1,NumOfComms);
for k = 1:NumOfComms
    PlotCommand{k} = varargin{2*k-1};
    PlotComVal{k} = varargin{2*k};
end %for k = 1:NumOfComms

%%% Setting Colors
TTotalBatsNumber = DataToAnalyze.AllParams.SimParams.TotalBatsNumber;
TotalPreysNumber = DataToAnalyze.AllParams.SimParams.TotalPreysNumber;

[colStruct] = myColorScheme(SimParams.TestMode, TotalBatsNumber, TotalPreysNumber);
ColorOfPreys     = colStruct.PreysPos;
ColorHuntingMan  = colStruct.Approach; % [ 0 1 0];  % Green
ColorObsMan      = colStruct.ObsMan; % [0 0.7 0]; % Green
ColorFindingPrey = colStruct.detections; %[1 0 1]; % magenta
ColorSearchMan = colStruct.ForagingMan; % [0 0 0]; %Black


ColorOfStages = {'k', 'r', 'b' ,[0.8 0.8 0],  [0.5 0.0 0.5]};



%%% Finding the Indices and plot

SampleTime = AllParams.SimParams.SampleTime;

MaxIndex = length(VectorToPlot);

for n =1:length(uStages)
    IndexSt.(uStages{n}) = getManuverInd(ReferenceCell, uStages(n), MaxIndex);
end % for n

ApproachIndex = getManuverInd(ReferenceCell, Value1, MaxIndex);
ObsManIndex = getManuverInd(ReferenceCell, Value2, MaxIndex);
SearchIndex = getManuverInd(ReferenceCell, Value3, MaxIndex);
BuzzIndex = getManuverInd(ReferenceCell, Value4, MaxIndex);


%%% The Plots

h2 = FigHandle;

if ~isempty(ApproachIndex)
    hp3 = plot(h2,ApproachIndex.*SampleTime, VectorToPlot(ApproachIndex),'o-', 'color', ColorHuntingMan ,'linewidth',2); % the position whenManuvering
% % %     hp3 = plot(h2,HuntingIndex.*SampleTime, VectorToPlot(HuntingIndex), 'color', 'm' ,'linewidth',2); % the position whenManuvering
    hold on
end % if ~isempty(HuntingIndex)

if ~isempty(ObsManIndex)
    hp1 = plot(h2,ObsManIndex.*SampleTime, VectorToPlot(ObsManIndex),'o-', 'color', ColorObsMan ,'linewidth',2); % the position whenManuvering
% % %     hp1 = plot(h2,ObsManIndex.*SampleTime, VectorToPlot(ObsManIndex), 'color', 'b' ,'linewidth',2); % the position whenManuvering
end % if ~isempty(HuntingIndex)

if ~isempty(SearchIndex)
    hp2 = plot(h2,SearchIndex.*SampleTime, VectorToPlot(SearchIndex),'o-', 'color', ColorSearchMan ,'linewidth',2); % the position whenManuvering
% % %      hp2 = plot(h2,ForagingIndex.*SampleTime, VectorToPlot(ForagingIndex), 'color', 'c' ,'linewidth',2); % the position whenManuvering
end % if ~isempty(HuntingIndex)
 
if ~isempty(BuzzIndex)
    hp4 = plot(h2,BuzzIndex.*SampleTime, VectorToPlot(BuzzIndex),'o-', 'color', ColorBuzzMan ,'linewidth',2); % the position whenManuvering
% % %      hp4 = plot(h2,ForagingIndex.*SampleTime, VectorToPlot(ForagingIndex), 'color', 'c' ,'linewidth',2); % the position whenManuvering
end % if ~isempty(HuntingIndex)
 

%%% Robust code for any stages 
% for nStage =1:length(uStages)
%     if ~isempty(IndexSt.(uStages{nStage}))
%         
%         for kPlotCom = 1:NumOfComms
%             hp(nStage) = plot( h2,...
%                 IndexSt.(uStages{nStage}) .* SampleTime, VectorToPlot(IndexSt.(uStages{nStage})),...
%                 'color', ColorOfStages{nStage} ,'linewidth',2);
%             set(hp(nStage), PlotCommand{kPlotCom} , PlotComVal{kPlotCom} );
%         end % for kPlotCom = 1:NumOfComms
% %         t = text(180,textYpos(kk), uStages{kk}, 'color', ColorOfStages{kk})
%     end % if ~isempty(IndexSt.(uStages(nStage)))
% %     IndexSt.(uStages(nStage)) = getManuverInd(ReferenceCell, uStages(kk), MaxIndex);
% end % for nStage
% L = legend((uStages{:}));    

%     for k = 1:NumOfComms
%         if ~isempty(ObsManIndex)
%             hp1 = plot(h2,ObsManIndex.*SampleTime, VectorToPlot(ObsManIndex), 'color', ColorObsMan ,'linewidth',2);
%             set(hp1, PlotCommand{k} , PlotComVal{k});
%         end
%         if ~isempty(ForagingIndex)
%             hp2 = plot(h2,ForagingIndex.*SampleTime, VectorToPlot(ForagingIndex), 'color', ColorForagingMan ,'linewidth',2);
%             set(hp2, PlotCommand{k} , PlotComVal{k});
%         end
%         if ~isempty(HuntingIndex)
%             hp3 = plot(h2,HuntingIndex.*SampleTime, VectorToPlot(HuntingIndex), 'color', ColorHuntingMan ,'linewidth',2);
%             set(hp3, PlotCommand{k} , PlotComVal{k});
%         end
%     end %for k = 1:NumOfComms

  
    
axis('auto');

end
%%%%%%%%%
% % %  movved to indepency


% % % function [IndVec] = getManuverInd(ManueverTypeCell,ManToSearch, MaxIndex)
% % % 
% % % %     IndexCell = strfind(ManueverTypeCell, ManToSearch);
% % % %     FullIndex = find(~(cellfun('isempty', IndexCell)));
% % %     IndexCell = strcmp(ManueverTypeCell, ManToSearch);
% % %     FullIndex = find(IndexCell);
% % %     
% % %     MaxVecIndex = find(FullIndex>=MaxIndex,1)-1;
% % %     
% % %     if ~isempty(MaxVecIndex)
% % %         IndVec = FullIndex(1:MaxVecIndex);
% % %     else % if ~isempty(MaxVecIndex)
% % %         IndVec = FullIndex;
% % %     end  % if ~isempty(MaxVecIndex)
% % % 
% % % 
% % % end % function [IndVec] = getManuverInd

