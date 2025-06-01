function [clusterCenters, numPointsInCluster, pointClusterIndices] = clusterPoints2(positions, maxDistThreshold)
    % Number of points
    numPoints = size(positions, 1);

    % Initialize cluster indices
    pointClusterIndices = zeros(numPoints, 1);
    currentCluster = 0;

    % Function to perform DFS and mark connected components
    function dfs(point, clusterIndex)
        stack = point;
        while ~isempty(stack)
            current = stack(end);
            stack(end) = [];
            if pointClusterIndices(current) == 0
                pointClusterIndices(current) = clusterIndex;
                distances = sqrt(sum((positions - positions(current, :)).^2, 2));
                neighbors = find(distances <= maxDistThreshold);
                stack = [stack; neighbors(pointClusterIndices(neighbors) == 0)];
            end
        end
    end

    % Identify clusters using DFS
    for i = 1:numPoints
        if pointClusterIndices(i) == 0
            currentCluster = currentCluster + 1;
            dfs(i, currentCluster);
        end
    end

    % Get the number of clusters
    numClusters = currentCluster;

    % Initialize output variables
    clusterCenters = zeros(numClusters, 2);
    numPointsInCluster = zeros(numClusters, 1);

    % Calculate the parameters for each cluster
    for i = 1:numClusters
        clusterMembers = positions(pointClusterIndices == i, :);
        clusterCenters(i, :) = mean(clusterMembers, 1);
        numPointsInCluster(i) = size(clusterMembers, 1);
    end
end