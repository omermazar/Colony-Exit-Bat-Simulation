function [numPointsInCluster, pointClusterIndices, clusterCenters] = clusterPoints(positions, minDistThreshold)
%%%% New June2024    
if size(positions,1) <= 1
    % no need for clustering
    numPointsInCluster = size(positions,1);
    pointClusterIndices = size(positions,1);
    clusterCenters = positions;
else
    % Perform hierarchical clustering
    Z = linkage(positions, 'single', 'euclidean');

    % Form clusters based on the distance threshold
    pointClusterIndices = cluster(Z, 'cutoff', minDistThreshold, 'criterion', 'distance');
    
    epsilon = minDistThreshold;
    minPts = 1;  % Minimum number of points to form a dense region
    pointClusterIndices = dbscan(positions, epsilon, minPts);
    
    % Get the number of clusters
    numClusters = max(pointClusterIndices);

    % Initialize output variables
    clusterCenters = zeros(numClusters, 2);
    numPointsInCluster = zeros(numClusters, 1);

    % Calculate the parameters for each cluster
    % c = lines(numClusters);
    if nargout >= 3
        for i = 1:numClusters
            clusterMembers = positions(pointClusterIndices == i, :);
            clusterCenters(i, :) = mean(clusterMembers, 1);
            numPointsInCluster(i) = size(clusterMembers, 1);
            % plot(positions(pointClusterIndices == i, 1), positions(pointClusterIndices == i, 2), 'o', ...
            %     'MarkerSize', 14, 'Color', c(i,:), 'MarkerFaceColor', c(i,:))
            % plot(clusterCenters(i, 1), clusterCenters(i, 2), '*', 'Color', 'k')
        end
    end % if nargout >= 3

end % if size(positions,1) > 1

% % % % Define a 2D matrix with some points
% % % positions = [
% % %     1, 2;
% % %     3, 4;
% % %     1.5, 2.5;
% % %     5, 6;
% % %     2, 2;
% % %     7, 8
% % % ];
% % % 
% % % % Define the minimum distance threshold
% % % minDistThreshold = 2.5;
% % % 
% % % % Cluster the points and get the desired parameters
% % % [clusterCenters, numPointsInCluster, pointClusterIndices] = clusterPoints(positions, minDistThreshold);
% % % 
% % % % Display the results
% % % disp('Cluster Centers (Average X-Y positions):');
% % % disp(clusterCenters);
% % % 
% % % disp('Number of points in each cluster:');
% % % disp(numPointsInCluster);
% % % 
% % % disp('Cluster indices for each point:');
% % % disp(pointClusterIndices);
