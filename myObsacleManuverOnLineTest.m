
function [ax]= myObsacleManuverOnLineTest(nStart, BAT, BatNum, CurrPulseNum, FoundObstacles, Terrain,ax, tetaCmd, maxDist, prev_xFinds, prev_yFinds, ...
    sortedXEst, sortedYEst, sortedAngles, sortedDistancesByAngles, ...
    edges, polynomes, wallIx, landMarks, exitPotentialsIx, selectedExitIx, wallOrMark, ManueverCmdStruct, currFoundObstacles)

% nTime = BAT(1).TransmittedPulsesStruct(AnalyzedPulseNum).StartPulseTime
% BatNum =1 
% CurrPulseNum = AnalyzedPulseNum
if nargin < 7
    fig = figure;
    ax = gca; hold on
else
%     ax = fig.Children(1);
end
if nargin <10
    prev_xFinds = [];
    prev_yFinds = [];
end 

%%% the colors
[colStruct] =  myColorScheme(); % (BatDATA.AllParams.SimParams.TestMode);

hold on
xyRes = 0.01;

%%% Bat Position
nEnd = min(nStart + BAT(BatNum).TransmittedPulsesStruct(CurrPulseNum).IPItoNextPulse-1, numel( BAT(BatNum).xBati));
xBat  = BAT(BatNum).xBati(nStart);
yBat  = BAT(BatNum).yBati(nStart);
tetaB = BAT(BatNum).Teta(nStart);
plot(ax, xBat*xyRes , yBat*xyRes, '.k', 'DisplayName','BatStart', 'MarkerSize', 9, 'LineWidth', 1)
% plot(ax, BAT(BatNum).xBati(nStart:nEnd)*xyRes , BAT(BatNum).yBati(nStart:nEnd)*xyRes, '.-k', 'DisplayName','manuever', 'MarkerSize', 9)
% quiver(xBat*0.01, yBat*0.01, 0.5*cos(tetaB),0.5*sin(tetaB), 'b', 'MaxHeadSize', 0.5, 'LineWidth', 1.5, 'DisplayName', 'movemnt')


%%%% Terrain
for k = 1:2:size(Terrain,1)-1
    plot(ax, Terrain(k,:)*xyRes,Terrain(k+1,:)*xyRes, '.-', 'Color', colStruct.Terrain, 'DisplayName', 'CaveWalls');
end
% plot(ax, BAT(BatNum).ObsInBeamStruct(CurrPulseNum).xOBSTC*xyRes, BAT(BatNum).ObsInBeamStruct(CurrPulseNum).yOBSTC*xyRes, 'o', 'Color', colStruct.Terrain, 'DisplayName','Pre-Detected Obs')

%%%% Detections
xDetected = (BAT(BatNum).xBati(nStart) + FoundObstacles.Distances.*cos(BAT(BatNum).Teta(nStart) + FoundObstacles.Angles) )*xyRes;
yDetected = (BAT(BatNum).yBati(nStart) + FoundObstacles.Distances.*sin(BAT(BatNum).Teta(nStart) + FoundObstacles.Angles) )*xyRes;
% plot(ax, xDetected, yDetected, 'o', 'LineWidth', 1.5, ...
%              'MarkerEdgeColor', colStruct.obsEstimate, 'MarkerSize', 6, 'DisplayName', 'Clustered Obs')
colStruct.obsEstimate = 'k'; %  % changed 18May25
plot(ax, xDetected, yDetected, 's', 'MarkerFaceColor', colStruct.obsEstimate, ...
                    'MarkerEdgeColor', colStruct.obsEstimate, 'MarkerSize', 6, 'DisplayName', 'Clustered Obs')

%%%% Memory
% plot(ax, prev_xFinds*xyRes, prev_yFinds*xyRes, '.', 'color', colStruct.obsEstimate, 'DisplayName', 'Memory detections', 'LineWidth', 1.5)
% witohout cluster
if isfield(BAT(BatNum), 'ObsInMemory')
    plot(BAT(BatNum).ObsInMemory(CurrPulseNum).xAll*xyRes, BAT(BatNum).ObsInMemory(CurrPulseNum).yAll*xyRes, ...
        '.',  'color', 'k', 'DisplayName', 'Memory all', 'LineWidth', 1.5, 'MarkerSize', 6)
end
%%%% only Current
plot(currFoundObstacles.xFindsEst*0.01, currFoundObstacles.yFindsEst*0.01, 'o', ...
    'color','g', 'DisplayName', 'Current', 'LineWidth', 1.5, 'MarkerSize', 9)
% plot(currFoundObstacles.xFindsEst*0.01, currFoundObstacles.yFindsEst*0.01, 'x', ...
%     'color', colStruct.obsEstimate, 'DisplayName', 'Current', 'LineWidth', 1.5, 'MarkerSize', 9)


if nargin >= 8
    tetaReq =  tetaCmd + BAT(BatNum).Teta(nStart);
    %%%% Add quiver for the command
    if ~isnan(tetaReq)
        % tetaCmd = BAT(BatNum).ManueverCmdStruct(CurrPulseNum+1).Angle2HuntedPrey;
        % tetaReq = BAT(BatNum).Teta(nStart) + tetaCmd;
        xStart = BAT(BatNum).xBatPos(nStart);
        yStart = BAT(BatNum).yBatPos(nStart);
        r = 0.5;
        quiver(xStart, yStart, r*cos(tetaReq), r*sin(tetaReq), 'MaxHeadSize', 1, ...
            'Color', 'k', 'DisplayName', 'tetaReq', 'LineWidth', 1.5)

        % plot(ax, 0.01*(BAT(BatNum).xBati(nEnd) + [0, maxDist* cos(directionCmd)]) , ...
        %     0.01*(BAT(BatNum).yBati(nEnd) + [0, maxDist* sin(directionCmd)]), ...
        %     '--m', 'DisplayName','Direction Commnand', 'MarkerSize', 6, 'LineWidth', 1.5)
    else 
        tempTeta = atan2(BAT(BatNum).yBati(nEnd) - BAT(BatNum).yBati(nStart), BAT(BatNum).xBati(nEnd) - BAT(BatNum).xBati(nStart));
        plot(ax, 0.01*(BAT(BatNum).xBati(nEnd) + [0, maxDist* cos(tempTeta)]) , ...
            0.01*(BAT(BatNum).yBati(nEnd) + [0, maxDist* sin(tempTeta)]), ...
            ':b', 'DisplayName','Direction Rand', 'MarkerSize', 6, 'LineWidth', 1.5)
        end % if isnan(directionCmd)

end % if nargin

%%%% add Decision guidelines
%%% Walls 

if ~isempty(edges)
    % THe selected wall
    pointsInWall = [edges(wallIx,1):edges(wallIx,2)];
    % The parameters of the wall
    m = polynomes(wallIx, 1); %slope
    b = polynomes(wallIx, 2); % intersecr
    x = sortedXEst(pointsInWall);
    y = sortedYEst(pointsInWall);
    % plot(x*0.01, y*0.01, 'b-s', 'LineWidth', 2, 'DisplayName', 'Wall', 'MarkerSize', 7)
    plot([x(1), x(end)]*0.01, (m*[x(1),x(end)]+b)*0.01, '--', 'color', 'b', ... % colStruct.obsEstimate, ...
        'LineWidth', 3, 'DisplayName', 'Wall', 'MarkerSize', 7)
    % plot(xw1*0.01, yw1*0.01, 'k*')
    % plot(xw2*0.01, yw2*0.01, 'r*')
    % 
    % all walls
    % for kWall = 1:size(edges,1)
    %     % dont plot te selected
    %     if kWall == wallIx
    %         continue
    %     end
    %     % else
    %     k1 = edges(kWall,1);
    %     k2 = edges(kWall,2);
    %     m = polynomes(kWall, 1); %slope
    %     b = polynomes(kWall, 2); % intersecr
    %     x = linspace(sortedXEst(k1), sortedXEst(k2),20);
    %     y = linspace(sortedYEst(k1), sortedYEst(k2),20);
    %     % plot(x*0.01, y*0.01, 'b-s', 'LineWidth', 1, 'DisplayName', 'Wall', 'MarkerSize', 4)
    %      plot([x(1), x(end)]*0.01, (m*[x(1),x(end)]+b)*0.01, '--', 'Color' , colStruct.obsEstimate, ...
    %          'LineWidth', 1, 'DisplayName', 'Wall', 'MarkerSize', 7)
    %     %  plot(sortedXEst([k1:k2])*0.01, m*sortedXEst([k1:k2])*0.01+ b*0.01, 'b-s', 'LineWidth', 1, 'DisplayName', 'Wall')
    % 
    %     % % dist1 = sqrt( (sortedXEst(k1)-x).^2+(sortedYEst(k1)-y).^2 );
    %     % % [d1, m1] = min(dist1);
    %     % % xx1 = x(m1); yy1=y(m1);
    %     % % 
    %     % % dist2 = sqrt( (sortedXEst(k2)-x).^2+(sortedYEst(k2)-y).^2 );
    %     % % [d2, m2] = min(dist2);
    %     % % xx2 = x(m2); yy2=y(m2);
    %     % % plot(x([m1,m2])*0.01, y([m1,m2])*0.01, 'b-s', 'LineWidth', 1, 'DisplayName', 'Wall')
    % end
end % if ~isempty(edges)

%%%%% landMarks
% if ~isempty(landMarks)
%     plot(landMarks(:,1)*0.01, landMarks(:,2)*0.01, '*g', 'MarkerSize', 20, 'LineWidth', 2, 'DisplayName', 'LandMarks')
%     if wallOrMark == "landMark"
%         xWall = sortedXEst(selectedExitIx(1));
%         yWall = sortedYEst(selectedExitIx(1));
%         xMark = landMarks(selectedExitIx(2),1);
%         yMark = landMarks(selectedExitIx(2),2);
%         plot([xWall, xMark]*0.01, [yWall, yMark]*0.01, 'og--', 'MarkerSize', 20, 'LineWidth', 2, 'DisplayName', 'Selected Gap')
%     end
% end
%%%%% caveExit potentials
exitEdges = selectedExitIx; % [selectedGap', exitPotentialsIx'+1];
% plot(sortedXEst(exitPotentialsIx)*0.01, sortedYEst(exitPotentialsIx)*0.01, 'db', 'MarkerSize', 12, 'DisplayName', 'Exit potential')
% plot(sortedXEst([selectedExitIx])*0.01, sortedYEst([selectedExitIx])*0.01, '*b--', ...
%     'LineWidth', 1.5, 'MarkerSize', 16, 'DisplayName', 'Suspected Exit')

%%%% BatAvoidance manuever - New June2024
if strcmp(ManueverCmdStruct.ManueverStage, 'AvoidBatMan')
    
    Dist2Bat  = ManueverCmdStruct.Dist2Bat;
    Angle2Bat = ManueverCmdStruct.Angle2Bat;
    xBatToAvoid = xBat + Dist2Bat*cos(Angle2Bat+tetaB);
    yBatToAvoid = yBat + Dist2Bat*sin(Angle2Bat+tetaB);

    plot([xBat, xBatToAvoid]*xyRes , [yBat, yBatToAvoid]*xyRes, 'o--r', 'MarkerSize', 12, 'LineWidth', 1.5, 'DisplayName', "batToAvoid")
end

% plot(ax, FoundObstacles.xFinds*xyRes, FoundObstacles.yFinds*xyRes, 'o', 'MarkerSize', 12, ...
%     'color', colStruct.Terrain, 'DisplayName', 'Detected Obs')
% legend( 'Location','bestoutside')
ax= gca;
plotPolar = 0;
if  plotPolar
    figure
    % all points
    polarplot(sortedAngles, sortedDistancesByAngles, 'og', 'MarkerSize', 6, 'DisplayName', "AllDetected")
    hold on
    
    % Walls
    if ~isempty(edges)
        for kWall = 1:size(edges,1)
                polarplot(sortedAngles(edges(kWall,1):edges(kWall,2)) , ...
                          sortedDistancesByAngles(edges(kWall,1):edges(kWall,2)), 'sb', 'MarkerSize', 4, 'DisplayName', "Wall")
        end
    end % edges

    % THe direction command
    polarplot([0,tetaCmd], [0, 200], 's-b', "DisplayName", "Required Direction")
    polarplot([0,0], [0, 200], 's-k', "DisplayName", "Flight direction")

    % the Exitgaps
    polarplot(sortedAngles(exitPotentialsIx), sortedDistancesByAngles(exitPotentialsIx), ...
         'db', 'MarkerSize', 12, 'DisplayName', 'Exit potential')
    polarplot(sortedAngles(selectedExitIx), sortedDistancesByAngles(selectedExitIx), ...
        '*b', 'MarkerSize', 12, 'DisplayName', 'Selected gap')

end
