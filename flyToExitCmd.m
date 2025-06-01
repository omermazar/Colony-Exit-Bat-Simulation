
function [ManueverCmdStruct, selectedGap, wallOrMark] = flyToExitCmd(ManueverCmdStruct, exitPotentialsIx, edgeToLandMarkPotential, BAT, AnalyzedPulseNum,  ...
                    sortedAngles, sortedDistancesByAngles, sortedXEst, sortedYEst, ...
                    landMarks, angleToLandMarks, distToLandMarks, ...
                    distToKeepFromWall, MinDistanceAllowed, DistanceToBuzz)

%%% New June2024
% flyToExitCmd
%%% Fly toward the middleof the closest gap to your direction
% the input angles are soterd from right to left

% [ManueverCmdStruct, selectedGap] = flyToExitCmd(ManueverCmdStruct, exitPotentialsIx, BAT, AnalyzedPulseNum,  ...
%                     sortedAngles, sortedDistancesByAngles, sortedXEst, sortedYEst, ...
%                     distToKeepFromWall, MinDistanceAllowed, DistanceToBuzz)

% returns the following varaibles:
% ManueverCmdStruct: the full Command to exit the cave (or not)
% selectedGap: the indices of corners of the selected gap (indices of sortedAngles)

% Find the closest Gap to the flight dierection
if any(exitPotentialsIx)
    edgesDierctions = sortedAngles(exitPotentialsIx);
    [~, minInd] = min(abs(edgesDierctions),[],"all");
    [row, col] = ind2sub(size(edgesDierctions), minInd);
    selectedGap = exitPotentialsIx(row,:);
    % fly to the middle of the gap
    requiredDirection = mean(sortedAngles(selectedGap));
    distToGap = min(sortedDistancesByAngles(selectedGap));
    wallOrMark = "Wall";
else
    % if no detections - ignore
    requiredDirection = inf;
    distToGap = inf;
    selectedGap = [];
    wallOrMark = "None";
end



% Find the Closest angle to the land marks
if any(edgeToLandMarkPotential)
    ixMark = edgeToLandMarkPotential(:,2);
    marksDirections = angleToLandMarks(ixMark);
    marksDists = distToLandMarks(ixMark);
    ixEdge  = edgeToLandMarkPotential(:,1);
    edgesMDierctions =  sortedAngles(ixEdge);
    [~, emInd] = min(abs(mean([edgesMDierctions', marksDirections],2)) );
    emGapDirection = mean([edgesMDierctions(emInd), marksDirections(emInd)]);
    emGapDist = mean([marksDists(emInd), sortedDistancesByAngles(ixEdge(emInd))]);
else
    emGapDist = inf;
    emGapDirection = inf;
end


%%% overirde wall edges with land-marks
% prefere walls
angleDiff = 10/180*pi;
distDiff = 50; %0.5 m
if (abs(emGapDirection) < pi/2-angleDiff ) && ...
        ( abs(emGapDirection) < abs(abs(requiredDirection)-angleDiff) ) && ...
    ( emGapDist < (distToGap - distDiff) )
    requiredDirection = emGapDirection;
    distToGap = emGapDist;
    selectedGap = edgeToLandMarkPotential(emInd,:);
    wallOrMark = "landMark";
end

%%% The Manuever Command
if ~isinf(requiredDirection)
    ManueverCmdStruct.ManueverStage    = 'CaveExit'; % NEW
    ManueverCmdStruct.ManueverType     = 'ExitMan'; % keep old definitions for minimize chnges in pther functions
    ManueverCmdStruct.Dist2HuntedPrey  = distToGap;
    ManueverCmdStruct.Angle2HuntedPrey = requiredDirection;
end

%%% Check if need to Buzz
if (ManueverCmdStruct.Dist2HuntedPrey <= DistanceToBuzz)
    ManueverCmdStruct.ManueverStage = 'Buzz';
end