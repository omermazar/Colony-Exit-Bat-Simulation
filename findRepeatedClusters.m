function [echoesPerCall, repeatedX, repeatedY, repeatedTargetInd, repeatedCenters, allX, allY] = findRepeatedClusters(Struct, echoesRatioTH, minDistThreshold)

%  %%%% New June2024
% retuturn the points which are part of clusters bigger than echoesTH/ them
% maximal distance in the clutters is defined by minDistThreshold;
% Inputs:
% Struct = BAT.ObsFindsStruct(minPulse:AnalyzedPulseNum-1)
% echoesRatioTH = 2;
% minDistThreshold = AllParams.TerrainParams.TerrainGrid*0.8/AllParams.SimParams.xyResolution; % th
% Outputs:
%   echoesPerCall - how many echoes returned from each call
%   [repeatedX, repeatedY] - the position od the detected cluster with more than the minimal reuired points per clusrter
%   repeatedTargetInd - the indices of the detections per cluster 
%   repeatedCenters =  [repeatedX, repeatedY]
%   [allX, allY]    = all cluster detected (including one points)


% extract the x-y positions of the detecteions

xIn = [Struct.xFindsEst]';
yIn = [Struct.yFindsEst]';
repeatedTargetInd = [Struct.DetectedPreyNum];
prevLocs = [xIn, yIn];

% find how many calls in the struct contain setections
echoesPerCall = cellfun(@numel, {Struct.xFindsEst});
numOfRelevantCalls = sum(echoesPerCall > 0);
if echoesRatioTH >= 1
    % take the PAramete fro the All PAram itself
    echoesTH = echoesRatioTH;
else
    % The parameter is the Ratio
    echoesTH = numOfRelevantCalls*echoesRatioTH;
end

if isempty(prevLocs)
    echoesPerCall = 0;
    repeatedX     = [];
    repeatedY     = [];
    repeatedTargetInd =[];
    repeatedCenters = [];

else
    %% Method 1 - full clustering with distance TH
    % clusters - a cell-array with all inpout points clustered in cells,
    % centers - XY matrix of the center of each cluster
    [clusters, centers, indicesCluster] = distanceTHClustering(prevLocs, minDistThreshold);
    numPointsInCluter = cellfun(@(x) size(x,1), clusters);
    % remove clusters with not enough points
    isRepeated = numPointsInCluter >= echoesTH;
    repeatedX         = centers(isRepeated,1)';
    repeatedY         = centers(isRepeated,2)';
    repeatedCenters   = centers(isRepeated,:)';
    repeatedTargetInd = indicesCluster(isRepeated);
    allX              = centers(:,1)';
    allY              = centers(:,2)';

    % %% Method 2 - filter by pairwise distcane (no clusters)
    % % matrix of pairise distances
    % dAll = pdist2(prevLocs, prevLocs);
    % % We sum the self-distance too, because the if the echoTH < 1 than one point is enough 
    % nClosePoints = sum(dAll < minDistThreshold,2); 
    % isRepeated = nClosePoints >= echoesTH ; % 
    % repeatedCenters = [];

    %%
    %%%% Output
    % % validLocs = prevLocs(isRepeated,:);
    % repeatedX =  prevLocs(isRepeated,1)';
    % repeatedY =  prevLocs(isRepeated,2)';
    % repeatedTargetInd = repeatedTargetInd(isRepeated);

    
end % if isempty(prevLocs)

