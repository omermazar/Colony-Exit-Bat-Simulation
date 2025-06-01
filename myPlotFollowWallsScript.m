figure
hold on
% bat pos
nTime = BAT.TransmittedPulsesStruct(AnalyzedPulseNum).StartPulseTime;
xBat  = BAT.xBati(nTime);
yBat  = BAT.yBati(nTime);
tetaB = BAT.Teta(nTime);


% all detected
plot(sortedXEst*0.01, sortedYEst*0.01, '.k', 'MarkerSize', 12)

% ExitPotential
exitEdges = [exitPotentialsIx', exitPotentialsIx'+1];
plot(sortedXEst(exitEdges)*0.01, sortedYEst(exitEdges)*0.01, 'og', 'MarkerSize', 24)
plot(sortedXEst([selectedExitIx,selectedExitIx+1])*0.01, sortedYEst([selectedExitIx,selectedExitIx+1])*0.01, '*g', 'MarkerSize', 12)

% Walls
plot(sortedXEst(pointsInWall)*0.01, m*sortedXEst(pointsInWall)*0.01+ b*0.01, 'b-s', 'LineWidth', 2)
plot(xw1*0.01, yw1*0.01, 'k*')
plot(xw2*0.01, yw2*0.01, 'r*')

% the movement direction
quiver(xBat*0.01, yBat*0.01, 0.5*cos(tetaB),0.5*sin(tetaB), 'b', 'MaxHeadSize', 0.5, 'DisplayName', 'movemnt')
% the dierction to the selected point
quiver(xBat*0.01, yBat*0.01, 0.5*cos(tetaB+sortedAngles(selectedExitIx)),0.5*sin(tetaB+sortedAngles(selectedExitIx)), 'r', 'MaxHeadSize', 0.2, 'DisplayName', 'cmd')
% the command dierction
quiver(xBat*0.01, yBat*0.01, 0.5*cos(tetaB+requiredDirection),0.5*sin(tetaB+requiredDirection), 'k', 'MaxHeadSize', 0.2, 'DisplayName', 'cmd')
plot(plannedPoint(1)*0.01, plannedPoint(2)*0.01, 'sb')

% the reference plot
xRef = sortedXEst(ixToCheck(minIx)); % sorry for the indices ...
yRef = sortedYEst(ixToCheck(minIx));
plot(xRef*0.01, yRef*0.01, 'sr')

% teheplanned point
plot(plannedPoint(1)*0.01, plannedPoint(2)*0.01, 'sg')
refDist= turnRadius;
plannedPoint = [xBat*0.01 + 0.01*refDist*cos(tetaB+requiredDirection), yBat*0.01 + 0.01*refDist*sin(tetaB+requiredDirection)];

% ObsAvoidance
plot(sortedXEst(ixInfornt)*0.01, sortedYEst(ixInfornt)*0.01, 'or', 'MarkerSize', 12)

% get Distance form fig
[x,y] = ginput; sqrt(sum([diff(x).^2, diff(y).^2],2))


figure
plot(sqrt(diff(sortedXEst).^2+diff(sortedYEst).^2),'o-')
hold on
plot(xlim, diffSortedDistTH*[1,1], 'r--')
plot(xlim, maxDistForExit*[1,1], 'r--')  

tic 
for ii = 1:10000
    [edges, polynomes] = findLines(sortedXEst, sortedYEst, mindistDiff, minNumOfPoints);
end
toc