function [] = MyTimeColorByManueverPlot(VectorToPlot, BAT, AllParams, FigHandle, varargin)

% plot  the Desired Vector in colors of the manuever
% inputs- VectorToPlot, ManueverTypeCell, SampleTime, FigHandle
%   VectorToPlot- the vector we plot
%   BAT = the dat Strucct from the fight
%   AllParmas- the stauct from the simulation - 
%   FigHandle- the handle of the figure
%   VectorType - 'Data' or 'Referenence'
%   optional - varagin - comman to line stlye' sould be in pairs:
%       example: [varargin(1)m,varargin(2)] - [Marker,'+'] 

NumOfExtraInputs = length(varargin);
NumOfComms = floor(NumOfExtraInputs/2);
PlotCommand = repmat({''},1,NumOfComms);
PlotComVal =  repmat({''},1,NumOfComms);
for k = 1:NumOfComms
    PlotCommand{k} = varargin{2*k-1};
    PlotComVal{k} = varargin{2*k};
end %for k = 1:NumOfComms

%%% Setting Colors
TotalBatsNumber = AllParams.SimParams.TotalBatsNumber;
kColor1 = linspace(0.5,1, TotalBatsNumber);
ColorOfBat = zeros(TotalBatsNumber,3);
ColorOfBat(:,3) = kColor1; % Blue

TotalPreysNumber = AllParams.SimParams.TotalPreysNumber;
kColor2 = linspace(0.5,1, TotalPreysNumber);
ColorOfPreys = zeros(TotalPreysNumber,3);
ColorOfPreys(:,3) = kColor2;
ColorOfPreys(:,1) = kColor2; % magenta

ColorHuntingMan = [ 0 1 0];  % Green
ColorObsMan = [0 0.7 0]; % Green
ColorFindingPrey = [1 0 1]; % magenta
ColorForagingMan = [1 1 0]; %Black

%%% Finding the Indices and plot
ManueverTypeCell = BAT.ManueverTypeCell;
SampleTime = AllParams.SimParams.SampleTime;

MaxIndex = length(VectorToPlot);

HuntingIndex = getManuverInd(ManueverTypeCell, 'Hunting', MaxIndex);
ObsManIndex = getManuverInd(ManueverTypeCell, 'ObsMan', MaxIndex);
ForagingIndex = getManuverInd(ManueverTypeCell, 'Foraging', MaxIndex);

% HuntingIndexCell = strfind(ManueverTypeCell,'Hunting');
% HuntingIndex = find(~(cellfun('isempty', HuntingIndexCell)));
% MaxHuntingIndex = find(HuntingIndex>=MaxIndex,1)-1;
% HuntingIndex = HuntingIndex(1:MaxHuntingIndex); 

% ObsManIndexCell = strfind(ManueverTypeCell,'ObsMan');
% ObsManIndex = find(~(cellfun('isempty', ObsManIndexCell)));
% ObsManIndex = ObsManIndex(1:MaxIndex);
% 
% ForagingIndexCell = strfind(ManueverTypeCell,'Foraging');
% ForagingIndex = find(~(cellfun('isempty', ForagingIndexCell)));
% ForagingIndex = ForagingIndex(1:MaxIndex);

% Plot Reference
% if VectorType  == 'Reference'
%     
%     RefernceVec = zeros(1, MaxTime/SampleTime+1);
%     RefernceVec(HuntingIndex) = 2;
%     RefernceVec(ObsManIndex) = 1;
%     
%     VectorToPlot = RefernceVec;
%     
% end % if VectorType  = 'Reference'    


%%% The Plots

h2 = FigHandle;

if ~isempty(HuntingIndex)
    hp3 = plot(h2,HuntingIndex .*SampleTime, VectorToPlot(HuntingIndex), 'color', ColorHuntingMan ,'linewidth',2); % the position whenManuvering
    hold on
end % if ~isempty(HuntingIndex)

if ~isempty(ObsManIndex)
    hp1 = plot(h2,ObsManIndex .*SampleTime, VectorToPlot(ObsManIndex), 'color', ColorObsMan ,'linewidth',2); % the position whenManuvering
end % if ~isempty(HuntingIndex)

if ~isempty(ForagingIndex)
    hp2 = plot(h2,ForagingIndex .*SampleTime, VectorToPlot(ForagingIndex), 'color', ColorForagingMan ,'linewidth',2); % the position whenManuvering
end % if ~isempty(HuntingIndex)

for k = 1:NumOfComms
    if ~isempty(ObsManIndex)
        try %set
            set(hp1, PlotCommand{k} , PlotComVal{k});
        catch % try %set
            endtry %set
            pop = 'set (hp1, PlotCommand{k} , PlotComVal{k}) '
        end% try
    end
    if ~isempty(ForagingIndex)
        set(hp2, PlotCommand{k} , PlotComVal{k});
    end
    if ~isempty(HuntingIndex)
        set(hp3, PlotCommand{k} , PlotComVal{k});
    end
end %for k = 1:NumOfComms

axis('auto');

end
%%%%%%%%%

function [IndVec] = getManuverInd(ManueverTypeCell,ManToSearch, MaxIndex)

    IndexCell = strfind(ManueverTypeCell, ManToSearch);
    FullIndex = find(~(cellfun('isempty', IndexCell)));
    MaxVecIndex = find(FullIndex>=MaxIndex,1)-1;
    if ~isempty(MaxVecIndex)
        IndVec = FullIndex(1:MaxVecIndex);
    else % if ~isempty(MaxVecIndex)
        IndVec = FullIndex;
    end  % if ~isempty(MaxVecIndex)


end % function [IndVec] = getManuverInd

