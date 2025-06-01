function [exitPotentialsInd, edgeToLandMarkPotential] = findExitByWalls(edges, polynomes, landMarks, ...
    sortedAngles, sortedDistancesByAngles, sortedXEst, sortedYEst, angleToLandMarks, distToLandMarks, ...
    maxDetectionDist, diffSortedDistTH, maxDistForExit, AllParams)

 % this function find exit potintion by the ditances between all eges of
 % detected walls. Edges is a nx2 matrix of the indices of edgedes of each
 % wall.
 
 mindistDiff = 0.45 / 0.01; % the minima distance between wall edeges
 minDistPoint = 0.15/ 0.01; % the minimal distace fro pbstace to a valis exit potential 
 maxPointsAllowed = 3; % how many point can be close to the gap
 vAngle = 15/180*pi;

 numWalls = size(edges,1);
 % get the XY postions of the Edges
 Xedges = sortedXEst(edges);
 Yedges = sortedYEst(edges);

 %delete all far adeges
 ixTooFAr = sortedDistancesByAngles(edges) >= maxDetectionDist;
 Xedges(ixTooFAr) = inf;
 Yedges(ixTooFAr) = inf;
 % edges(ixTooFAr) = inf;
 
 % INIT
 exitPotentialsInd = [];
 edgeToLandMarkPotential= [];

 % for each row, Chek the distance condition with closest edge of other walls
 % in each row the reference point is the one closest (in Angle) to the
 % flight direction
 for k = 1:(numWalls-1)
     [~, refPointInd] = min(abs(sortedAngles(edges(k,:))) );
     refAngle = sortedAngles(edges(k,refPointInd));
     x1 = sortedXEst(edges(k,:));
     y1 = sortedYEst(edges(k,:));
     refX = Xedges(k, refPointInd);
     refY = Yedges(k, refPointInd);
     for j = (k+1):numWalls
         % distances = sqrt( (Xedges(k,:)'- Xedges(j,:)).^2 + (Yedges(k,:)'- Yedges(j,:)).^2);
         % [minDist, ii] = min(distances, [], "all");
         % if (minDist >= diffSortedDistTH) && (minDist <= maxDistForExit)
         %     [row, col] = ind2sub(size(distances), ii);
         %     exitPotentialsInd = [exitPotentialsInd; [edges(k,row), edges(j,col)] ];
         % end
         % if the walls areto close than no gap
         x2 = sortedXEst(edges(j,1):edges(j,2)); % sortedXEst(edges(j,:));
         y2 = sortedYEst(edges(j,1):edges(j,2)); % sortedYEst(edges(j,:));
         if any(pdist2([x1',y1'], [x2',y2']) < mindistDiff, "all")
             continue
         end

         % Distance to the nearest angle to the egde  in the other wall
         [~, curPiontInd] = min( abs( sortedAngles(edges(j,:))- refAngle) );
         curX = Xedges(j, curPiontInd);
         curY = Yedges(j, curPiontInd);
         minDist = sqrt( (curX-refX).^2 + (curY'- refY).^2);
         if (minDist >= diffSortedDistTH) && (minDist <= maxDistForExit)
             curEdge = [edges(k,refPointInd), edges(j,curPiontInd)];
             exitPotentialsInd = [exitPotentialsInd; curEdge ];
         end
     end % for j

     %%  %% Look for exit potials betwee edges and landMarks
     if any(landMarks)
         %  the distnce betwen the edge and the lansmarks
         distToLandMarks = sqrt((landMarks(:,1)- refX).^2 + (landMarks(:,2)- refY).^2);
         indValidDist =  find((distToLandMarks >= diffSortedDistTH) & (distToLandMarks <= maxDistForExit));
         % select the landMark that is closer to the ages
         if any(indValidDist)
             [~, mInd] = min(distToLandMarks(indValidDist));
             newPot = [edges(k,refPointInd), indValidDist(mInd)];
             edgeToLandMarkPotential = [edgeToLandMarkPotential; newPot];
         end

     end % if any(landMarks)
 end % for k

 %% select all points that are not edges
 ixNotEdge = ~ismember(1:numel(sortedXEst), reshape(edges, 1, []));
 xObs = sortedXEst(ixNotEdge);
 yObs = sortedYEst(ixNotEdge);

%% remove edgeToLandMarkPotential that are too close other obstacls
if ~isempty(edgeToLandMarkPotential)
    numPots  = size(edgeToLandMarkPotential,1);
  
    % init with all valid
    validPot = true(1,numPots);

    for ke = 1:numPots
        kStart = edgeToLandMarkPotential(ke,1);
        kEnd   = edgeToLandMarkPotential(ke,2);
        xExit = [sortedXEst(kStart), landMarks(kEnd,1)];
        yExit = [sortedYEst(kStart), landMarks(kEnd,2)];
        xlineExit = linspace(xExit(1), xExit(2), 10);
        ylineExit = linspace(yExit(1), yExit(2), 10);
        % remove the points that form the landMark
        distObsToLandMark = sqrt((xObs- xExit(2)).^2 + (yObs- yExit(2)).^2);
        ixCurr = distObsToLandMark <= 15;
        
        % Check wetherer there are obstacle detections between the current exitPotential
        distAllToLine =  pdist2([xObs(~ixCurr)',yObs(~ixCurr)'], [xlineExit',ylineExit']);
        distToLine = min(distAllToLine, [], 2);
        numCloseObsToLine = sum(distToLine < minDistPoint); % at least one point is 
        if numCloseObsToLine >= maxPointsAllowed-1
            validPot(ke) = false;
            continue; % one close wall is enough
        end

        %%% remove land Marks after detectition points of points and Before
        %%% Walls (Corners)
        isBetweenPoints = ( sortedAngles >= wrapToPi(angleToLandMarks(kEnd) - vAngle) ) & ...
                          ( (sortedAngles <= wrapToPi(angleToLandMarks(kEnd) + vAngle) ));
        isInfront = isBetweenPoints & (sortedDistancesByAngles < distToLandMarks(kEnd)); 
        if sum(isInfront) >= 2
            validPot(ke) = false;
        end
        
        %%% remove exitPotential that are infront   walls
        minAngle = min([sortedAngles(kStart), angleToLandMarks(kEnd)]);
        maxAngle = max([sortedAngles(kStart), angleToLandMarks(kEnd)]);
        distToExit = min([sortedDistancesByAngles(kStart), distToLandMarks(kEnd)]);
        isBetweenPoints = ( sortedAngles > minAngle) & (sortedAngles < maxAngle ); 
        isInfront = isBetweenPoints & (sortedDistancesByAngles < distToExit);
        isBehind  = isBetweenPoints & (sortedDistancesByAngles > distToExit) & (sortedDistancesByAngles < distToExit +100);
        if sum(isInfront) >=3 || sum(isBehind) >= 5
            validPot(ke) = false;
        end    
    end % for ke

    % update
    edgeToLandMarkPotential = edgeToLandMarkPotential(validPot,:);
end % if ~isempty(edgeToLandMarkPotential)



%% if there are not exit potentian - end
% dealing with 'broken walls' due to seperation
if AllParams.BatSonarParams.ObsClusterFlag
    fixWallsByClosedetections = true;
else
    fixWallsByClosedetections = false;
end

fixSplitedWallsFlag = true;


if isempty(exitPotentialsInd) || ~(fixSplitedWallsFlag || fixWallsByClosedetections) 
    return
end

%% Ignore Exit Potentials behind the Bat
angleTH = 100/180*pi;
ixBehind =  abs(sortedAngles(exitPotentialsInd(:,1))) > angleTH | abs(sortedAngles(exitPotentialsInd(:,2))) > angleTH;
exitPotentialsInd = exitPotentialsInd(~ixBehind,:);

%% Remove exitPotentialsInd that are too close to other walls or obstalce

%
numPots  = size(exitPotentialsInd,1);
% all detections thar are not edges
ixNotEdge = ~ismember(1:numel(sortedXEst), reshape(edges, 1, []));
xObs = sortedXEst(ixNotEdge);
yObs = sortedYEst(ixNotEdge);

% init wit all valid
validPot = true(1,numPots);

for ke = 1:numPots
    kStart = exitPotentialsInd(ke,1);
    kEnd   = exitPotentialsInd(ke,2);
    xExit = sortedXEst([kStart, kEnd]);
    yExit = sortedYEst([kStart, kEnd]);
    xlineExit = linspace(xExit(1), xExit(2), 10);
    ylineExit = linspace(yExit(1), yExit(2), 10);
    
    % Check wetherer there are obstacle detections between the current exitPotential
    if fixWallsByClosedetections
        distAllToLine =  pdist2([xObs',yObs'], [xlineExit',ylineExit']);
        distToLine = min(distAllToLine, [], 2);
        numCloseObsToLine = sum(distToLine < minDistPoint);
        if numCloseObsToLine >= maxPointsAllowed
            validPot(ke) = false;
            continue; % one close wall is enough
        end
    end % if fixWallsByClosedetections
    
    % Check the distances from this line to all other walls
    if fixSplitedWallsFlag
        for jWall = 1:numWalls
            % ignore walls thar are part of the currnet exit
            if any(ismember(edges(jWall,:), [kStart, kEnd]), "all")
                continue
            end
            % the wall cornnes
            x2 = sortedXEst(edges(jWall,:));
            y2 = sortedYEst(edges(jWall,:));
            % the distances between the wakkConers and Exit line
            lineDists =  pdist2([x2',y2'], [xlineExit',ylineExit']);
            % if the wall is closer than expected - delete the exit potential
            if any(lineDists < mindistDiff,"all")
                validPot(ke) = false;
                continue; % one close wall is enough
            end %
        end % for jWall
    end  % if fixSplitedWallsFlag

end %for ke
% the deletion
exitPotentialsInd = exitPotentialsInd(validPot,:);

%% for edges thar are common to several exit potenrtial select the exit with shorter gap (by angle)
if size(exitPotentialsInd,1) > 1
    m = size(exitPotentialsInd,1);  % Get the size of the matrix
    unique_indices = [];  % To store unique indices
    rows_to_keep = true(m, 1);  % Logical array to mark rows to keep

    % Loop through each row to check for duplicate indices
    for i = 1:m
        current_indices = exitPotentialsInd(i, :);

        % Check for duplicates
        duplicate_rows = find(any(exitPotentialsInd == current_indices(1), 2) | any(exitPotentialsInd == current_indices(2), 2));

        if length(duplicate_rows) > 1
            % Calculate the absolute difference for each duplicate row
            differences = sum(abs(exitPotentialsInd(duplicate_rows, 1) - exitPotentialsInd(duplicate_rows, 2)), 2);

            % Find the index of the row with the minimum difference
            [~, min_index] = min(differences);

            % Mark rows except the one with the minimum difference for removal
            rows_to_keep(duplicate_rows) = false;
            rows_to_keep(duplicate_rows(min_index)) = true;
        end
    end

    % keep only unique indices with minimum difference
    exitPotentialsInd = exitPotentialsInd(rows_to_keep, :);

end %  if size()

%% Remove exit Potential that are almost in the same line from the bat if they are not infornt of it
minAngleDiff = 2/180*pi; % 2deg
angleTH = 45/180*pi;
ixCloseAngles = abs(sortedAngles(exitPotentialsInd(:,1))- sortedAngles(exitPotentialsInd(:,2))) < minAngleDiff;
ixBehind =  abs(sortedAngles(exitPotentialsInd(:,1))) > angleTH |  abs(sortedAngles(exitPotentialsInd(:,2))) > angleTH;
ixRemove = ixCloseAngles & ixBehind;
exitPotentialsInd =exitPotentialsInd(~ixRemove,:);

%% remove ExitPotentials that are behind other obstacles
if ~isempty(exitPotentialsInd)
    validPot = true(size(exitPotentialsInd,1),1);
    for k =1:size(exitPotentialsInd,1)
        minAngle = min(sortedAngles(exitPotentialsInd(k,:)));
        maxAngle = max(sortedAngles(exitPotentialsInd(k,:)));
        distToExit = mean(sortedDistancesByAngles(exitPotentialsInd(k,:)));
        isBetweenPoints = ( sortedAngles > minAngle) & (sortedAngles < maxAngle ) & (sortedDistancesByAngles < distToExit);
         if sum(isBetweenPoints) >= 3
            validPot(k) = false;
        end
    end
    exitPotentialsInd = exitPotentialsInd(validPot,:);

end


 % % relEdges = allEdges(ixRel); 

 % % reshape all edged to one row
 % % allEdges = reshape(edges,1,[]);
 % % 
 % % % delete all far adeges
 % % ixRel = sortedDistancesByAngles(allEdges) <= maxDetectionDist;
 % % relEdges = allEdges(ixRel); 
 % % numPoints = sum(ixRel);
 % 
 % % the position of the relevant edges
 % Xedges = sortedXEst(relEdges);
 % Yedges = sortedYEst(relEdges);
 % 
 % % the distances between the relevnt wall-eges
 % distances = sqrt( (Xedges'- Xedges).^2 + (Yedges'- Yedges).^2);
 % % the diagonal to inf
 % distances(1:(numPoints+1):end) = Inf;
 % % take only upper triangular
 % distances = triu(distances);
 % 
 % % Check Distance condition
 % [rows, cols] = find(distances >= diffSortedDistTH & distances <= maxDistForExit);
 % detEges = [allEdges(rows)', allEdges(cols)'];
 % 
 % % remove edges of the same-walls
 % exitPotentialsInd = ~ismember(detEges, edges,"rows");
 % 
 