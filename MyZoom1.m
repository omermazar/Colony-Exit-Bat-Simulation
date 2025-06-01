
function[] = MyZoom1(obj, event_obj,HxyPlot, HTimePlot, DataToPlot, DataToAnalyze,...
                        SimParams, TerrainParams, Terrain)
% ----ZOOM LINKING Between figures
    ZoomFlag = 1;
    CurrentAx = event_obj.Axes;
    xLimits = xlim(CurrentAx);
    yLimits = ylim(CurrentAx);
    hLine = findobj(CurrentAx.Children,'Type','line');
% % %     hLine = handles.xyPlot.handleXYBatPosLine;
    hBar = findobj(CurrentAx.Children,'Type','bar');
    if ~isempty(hLine)    
        XD = hLine(1).XData;
        YD = hLine(1).YData;
    elseif ~isempty(hBar)  % ~isempty(h1)
        XD = hBar(1).XData;
        YD = hBar(1).YData;
    else % ~isempty(h1) % ~isempty(h2)
        XD =[];
        YD=[];
    end % ~isempty(h1)
%     
%     XD = CurrentAx.Children.XData;
%     YD = CurrentAx.Children.YData;
    IDnXYPlot = find(XD < xLimits(2) & XD > xLimits(1)...
        & YD < yLimits(2) & YD > yLimits(1));
    IDnTPlot = find(XD < xLimits(2) & XD > xLimits(1));
    if CurrentAx == HxyPlot
        MyTimePlot(HTimePlot,DataToAnalyze,DataToPlot,...
            SimParams, TerrainParams, Terrain,ZoomFlag, IDnXYPlot);
    elseif  CurrentAx == HTimePlot % if CurrentAx 
       MyXYPlot(HxyPlot,DataToAnalyze,DataToPlot,...
            SimParams, TerrainParams, Terrain, ZoomFlag, IDnTPlot);
    end % if CurrentAx 

