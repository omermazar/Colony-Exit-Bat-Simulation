function [landMarks, pointsPerMark] = findLandMarks(memFindsStruct, minDetectionMark, maxDistanceMark, edges, sortedXEst, sortedYEst) 

%  [landMarks] = findLandMarks(memFindsStruct, edges, minDetectionMark, maxDistanceMark);
%
% returns clouds of points with at least minDetectionMark detections and
% maximau distance of maxDistanceMar/ Lad Marks are not port of Walls


%% Input
detStruct.xFindsEst = [memFindsStruct.xFindsEst];
detStruct.yFindsEst = [memFindsStruct.yFindsEst];
detStruct.DetectedPreyNum = [memFindsStruct.DetectedPreyNum];

if ~exist('edges','var')
    edges = [];
end
% the minimum distance from land Mark to the line
minDistToLine = maxDistanceMark;

%% Remove detections realted to walls
if any(edges)
    % init
    nWalls = size(edges,1);
    indCloseDetection = [];
    % find valid point - far enough from the walls
    for k=1:nWalls
        ixWall = edges(k,1):edges(k,2);
        xWall = sortedXEst(ixWall);
        yWall = sortedYEst(ixWall);
        
        distAllToLine =  pdist2([detStruct.xFindsEst',detStruct.yFindsEst'], [xWall',yWall']);
        distToLine = min(distAllToLine, [], 2);
        indClose = find(distToLine < minDistToLine);
        indCloseDetection= unique([indCloseDetection; indClose]);
    end % for k
    % update the detection struct
    detStruct.xFindsEst(indCloseDetection)       = [];
    detStruct.yFindsEst(indCloseDetection)       = [];
    detStruct.DetectedPreyNum(indCloseDetection) = [];
else
    % dont remove
end %  if any(edges)

%% find the landMarks
if any(detStruct.xFindsEst)
    [~, xFinds, yFinds, repeatedInd] = findRepeatedClusters(detStruct, minDetectionMark, maxDistanceMark);
    landMarks = [xFinds', yFinds'];
    pointsPerMark = cellfun(@numel, repeatedInd);
else
    landMarks = [];
    pointsPerMark = [];
end

