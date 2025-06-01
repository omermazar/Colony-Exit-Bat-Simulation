function [clusters, centers, indCell] = distanceTHClustering(X, maxDist)
    % indCell is a cell with the indices of the input data in each clusters
    
    % Initialize clusters
    clusters = {};
    centers  = [];
    indCell = {}; 
    numPoints = size(X, 1);
    usedPoints = false(numPoints, 1);
    
    % Perform clustering
    while any(~usedPoints)
        % Find the first unused point
        idx = find(~usedPoints, 1);
        
        % Initialize a new cluster with the current point
        currentPoint = X(idx, :);
        cluster = currentPoint;
        usedPoints(idx) = true;
        
        % Create a list of points to check
        pointsToCheck = X(~usedPoints, :);
        distances = sqrt(sum((pointsToCheck - currentPoint) .^ 2, 2));
        
        % Find points within the maxDist from the current point
        neighbors = distances <= maxDist;
        
        % Add neighbors to the cluster
        cluster = [cluster; pointsToCheck(neighbors, :)]; %#ok<AGROW>
        clusterCenter = mean(cluster, 1);
        % Mark these neighbors as used
        unusedIdx = find(~usedPoints);
        usedPoints(unusedIdx(neighbors)) = true;
        
        % Add the new cluster to the list
        clusters{end+1} = cluster; %#ok<AGROW>
        centers = [centers; clusterCenter];
        indCell{end+1} = vertcat(idx, unusedIdx(neighbors));
    end
end

% Example usage
% % % Generate example data
% % % numPoints = 50;
% % % X = [300 + 200*rand(numPoints, 1), 50 + 250*rand(numPoints, 1)];
% % % 
% % % Set maximum distance from the cluster center
% % % maxDist = 10.0;
% % % 
% % % Perform recursive distance-based clustering
% % % clusters = efficientRecursiveDistanceBasedClustering(X, maxDist);
% % % 
% % % Plot the clusters
% % % figure;
% % % hold on;
% % % colors = lines(length(clusters));
% % % for i = 1:length(clusters)
% % %     scatter(clusters{i}(:,1), clusters{i}(:,2), 36, 'filled', 'MarkerFaceColor', colors(i,:));
% % % end
% % % title('Efficient Recursive Distance-Based Clustering');
% % % xlabel('X');
% % % ylabel('Y');
% % % hold off;
