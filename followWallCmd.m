function [ManueverCmdStruct, edges, polynomes, wallIx] = followWallCmd(ManueverCmdStruct, BAT, AnalyzedPulseNum,  ...
    sortedAngles, sortedDistancesByAngles, sortedXEst, sortedYEst, ...
    distToKeepFromWall, MinDistanceAllowed, DistanceToBuzz)

%%% New June2024
% Follow Wall command
%%% React to the obstale - Fly along the obstacle and Searh for exit
%%%% Turn to the direction with the further direction and keep
%%%% distance constant from wall
% the input angles are soterd from right to left

% [ManueverCmdStruct, edges, polynomes, wallIx] = followWallCmd(ManueverCmdStruct, BAT, AnalyzedPulseNum,  ...
%     sortedAngles, sortedDistancesByAngles, sortedXEst, sortedYEst, ...
%     distToKeepFromWall, MinDistanceAllowed, DistanceToBuzz) 
% returns the foolwing varaibles:
% ManueverCmdStruct: the full Command to foolw the wall (or not)
% edges: the indices of edges of the detected walls (indices of sortedAngles)
% polynomes: the parameteres of the linear fit for each waam (m,b, y = m*x+b)
% wallIx = the index of the wall that the bat selected to follow

% Parameters to follow wall
% distToKeepFromWall = 0.8;
minAngleDiff = 22.5 / 180 *pi; %( former - exitClearAngles)
mindistDiff = 0.45 / 0.01; % 0.4 % maximum distace betwen point to be condisered as oe wall (should be a leatts 2*grid)
inFrontAngleLimit = pi/4; % 45 dergrees
minNumOfPoints = 3;% 2; % the minimal number of points in a 'wall'
turnRadius = 1.5 / 0.01;

% The curent position of the bat
nTime = BAT.TransmittedPulsesStruct(AnalyzedPulseNum).StartPulseTime;
xBat  = BAT.xBati(nTime);
yBat  = BAT.yBati(nTime);
tetaB = BAT.Teta(nTime);

%%% In the navigation, consider only objects infront the bat
% ixValid = abs(sortedAngles) < pi/2;
% sortedAngles = sortedAngles(ixValid);
% sortedDistancesByAngles = sortedDistancesByAngles(ixValid);
% sortedXEst = sortedXEst(ixValid);
% sortedYEst = sortedYEst(ixValid);

% If less than three points- no wall
if numel(sortedAngles) <  minNumOfPoints
    edges    = [];
    polynomes = []; 
    wallIx   = [];

    return
end
%% cluster to walls and calculate directions of walls 
% follow walls by left and right angles AND by gaps in sorted angles 
% gaps may inidicate other walls
[edges, polynomes] = findLines(sortedXEst, sortedYEst, mindistDiff, minNumOfPoints);

%% split too long wall and check for different directions in walls
tooManyPointsTH = 10;
numPointsInWall = diff(edges,[],2);
indToCheck = find(numPointsInWall > tooManyPointsTH);
if ~isempty(indToCheck)
    tempEdges = edges; % temp variable
    tempPolys = polynomes;
    % follow the numbers of ne edges
    shiftInd = 0;
    for nCheck = indToCheck(:)'
        xWall = sortedXEst(edges(nCheck,1):edges(nCheck,2));
        yWall = sortedYEst(edges(nCheck,1):edges(nCheck,2));
        % refline = polynomes(nCheck);
        [newEdges, newlines, newPolys] = splitWallLines(xWall, yWall, minNumOfPoints, false);
        % replace current Wall with new Edges an polynomes 
        if size(newEdges,1) > 1
            shiftedEdges = edges(nCheck,1) + newEdges - 1;
            indToReplace = nCheck + shiftInd;
            tempEdges = [tempEdges((1:indToReplace-1),:); shiftedEdges; tempEdges((indToReplace+1:end),:)];
            tempPolys = [tempPolys((1:indToReplace-1),:); newPolys; tempPolys((indToReplace+1:end),:)];
            shiftInd = shiftInd + size(newEdges,1) -1;
        end
    end % for nCheck
    % replace with the new
    edges = tempEdges;
    polynomes = tempPolys;

end% if

%% TBD

%% Mnuever Command
% edges = [];
% if there are no estimated walls - avoid the nearest angle point
if isempty(edges)
    %% no walls
    % Check whether to fly  right or left to the reference point (+1) is counter-clockwise
    % if the the referenecre point is the 'left' - (+1)
    % if the the referenecre point is the 'right' -(-1)
    wallIx = [];
    polynomes = [];
    ixToCheck = [1, numel(sortedAngles)];
    anglesToCheck = sortedAngles(ixToCheck);
    % follow the closet angle-point to the bat's flight direction
    [~, minIx] = min(abs(anglesToCheck));
    refAngle = anglesToCheck(minIx);
    refDist = sortedDistancesByAngles(ixToCheck(minIx)); % sorry for the indices ...
    
    if minIx == 1 % most right
        signCorrect = -1;
    elseif minIx == numel(anglesToCheck) % most left
        signCorrect = +1;
    else
    end

    % follow the wall in cont distance
    requiredDirection = atan((refDist*sin(refAngle)+signCorrect*distToKeepFromWall) ./ (refDist*cos(refAngle)) );
    % the same sign as the original angle
    requiredDirection = sign(refAngle)*requiredDirection;
    plannedPoint = [xBat + refDist*cos(tetaB+requiredDirection), yBat + refDist*sin(tetaB+requiredDirection)];
    followObsFlag = true;

else
    %% Walls 
    % follow the wall with the minimal view angle infront of the bat
    ixToCheck = unique(edges(:));
    anglesToCheck = sortedAngles(ixToCheck);
    [~, ind0] = min(abs(anglesToCheck));
    wallIx = find(any(edges == ixToCheck(ind0),2));
    pointsInWall = [edges(wallIx,1):edges(wallIx,2)];

    % The parameters of the wall
    m = polynomes(wallIx, 1); %slope
    b = polynomes(wallIx, 2); % intersecr
    % Find the edgepoints of the wall end estimate its direction
    xw1 = sortedXEst(pointsInWall(1));  
    xw2 = sortedXEst(pointsInWall(end));
    yw1 = sortedYEst(pointsInWall(1)); %yw1 = m*xw1+b;
    yw2 = sortedYEst(pointsInWall(end)); %yw2 = m*xw2+b; 
    tetaW1 = sortedAngles(pointsInWall(1));
    tetaW2 = sortedAngles(pointsInWall(end));

    %%% find the closest point on the estimated wall to the edges - it's
    % %%% somtimes fall becase the vertical lines of walls 
    % % the line
    % x = linspace(sortedXEst(pointsInWall(1)), sortedXEst(pointsInWall(end)),20);
    % y = m*x+b;
    % dist1 = sqrt( (xw1-x).^2+(yw1-y).^2 );
    % [d1, m1] = min(dist1);
    % xx1 = x(m1); yy1=y(m1);
    % 
    % dist2 = sqrt( (xw2-x).^2+(yw2-y).^2 );
    % [d2, m2] = min(dist2);
    % xx2 = x(m2); yy2=y(m2);
    
    xx1 = xw1; yy1 = yw1;
    xx2 = xw2; yy2 = yw2;

    % wallDirection = atan2(yw2-yw1, xw2-xw1);
    % THe bat flies to the farther point in the wall
    % for long memory -  we should ignore far away walls and behind
    % the disance - 
    dw1 = sqrt((yBat-yy1)^2++(xBat-xx1)^2);
    dw2 = sqrt((yBat-yy2)^2++(xBat-xx2)^2);
    
    %%% THe decision
    % first go the edge infront and not the edge bhihind
    if (abs(tetaW1) < pi/2) && (abs(tetaW2) >= pi/2)
        flyTowardCmd = 1;
    elseif (abs(tetaW2) < pi/2) && (abs(tetaW1) >= pi/2)
        flyTowardCmd = 2;
    % if both are infront or bhenid - fly to the farther
    elseif dw1 <= 0.9*dw2 % fly twoard point 2
        flyTowardCmd = 2; 
       
    elseif dw2 <= 0.9*dw1 % fly toward point1
        flyTowardCmd = 1;
    
    % if distancesa  are alike - fly to the easier
    else % follow by closer angle

        if abs(tetaW1) < abs(tetaW2) % fly towardpoint 1
            flyTowardCmd = 1; 
        else % fly toward pont 2
            flyTowardCmd = 2; 
        end
    end
    
    if flyTowardCmd == 2
        wallDirection = atan2(yw2-yw1, xw2-xw1);
        xPlanned = xx2 + distToKeepFromWall*cos(pi/2 +wallDirection);
        yPlanned = yy2 + distToKeepFromWall*sin(pi/2 +wallDirection);
    else
        wallDirection = atan2(yw1-yw2, xw1-xw2);
        xPlanned = xx1 + distToKeepFromWall*cos(3*pi/2 +wallDirection);
        yPlanned = yy1 + distToKeepFromWall*sin(3*pi/2 +wallDirection);
    end
    
    % requiredDirection =  wrapToPi(wallDirection - tetaB);
    plannedPoint = [xPlanned, yPlanned];
    requiredDirection = wrapToPi(atan2(yPlanned-yBat, xPlanned-xBat)-tetaB);
    followObsFlag = true;
    

    refDist = turnRadius;
    % plannedPoint = [xBat + refDist*cos(tetaB+requiredDirection), yBat + refDist*sin(tetaB+requiredDirection)];
    

end % if ~isempty(edges)

% implented in  a different function
%% Look for the closest obstable infront of the bat and adjust flight path - Romoved to different function
% ixInfornt = abs(sortedAngles) <= inFrontAngleLimit;
% % the closest distance to the planned flight path
% x0 = sortedXEst(ixInfornt);
% y0 = sortedYEst(ixInfornt);
% ang0 = sortedAngles(ixInfornt);
% %%% caulcate distances between all point to the planned path
% % the palnned path
% x1 = xBat; y1= yBat;
% x2 = plannedPoint(1); y2 = plannedPoint(2);
% xx = linspace(x1,x2,10);
% yy = linspace(y1,y2,10);
% allDists = min(pdist2([x0',y0'],[xx',yy']),[],2);
% [minDist, mi] = min(allDists);
% % Calculate the distance using the formula
% % numerator = abs((y2-y1).*x0 - (x2-x1).*y0 + x2.*y1 - y2.*x1);
% % denominator = sqrt((y2-y1).^2 + (x2-x1).^2);
% % d = numerator ./ denominator;
% % minDist   = min(d);
% % if too close slow down and correct the angle
% if minDist < distToKeepFromWall
%     ManueverCmdStruct.ManAccelCommand = 'SlowDown';
% 
%     requiredDirection = requiredDirection - sign(ang0(mi))*pi/15;
% end
% 
% % check for minimal Distance - Crush
% if minDist < 2*MinDistanceAllowed
%     followObsFlag = false;
% else
%     % followObsFlag = true;
% end

%% the output
if followObsFlag
    ManueverCmdStruct.ManueverStage   = 'Approach'; %
    ManueverCmdStruct.ManueverType    = 'FollowWall'; % NEW
    ManueverCmdStruct.Dist2HuntedPrey = refDist;
    ManueverCmdStruct.Angle2HuntedPrey = requiredDirection;
    % Check if need to Buzz
    if (ManueverCmdStruct.Dist2HuntedPrey <= DistanceToBuzz)
        ManueverCmdStruct.ManueverStage = 'Buzz';
    end
end %if followObsFlag
