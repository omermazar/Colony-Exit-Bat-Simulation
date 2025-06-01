function [clusters, C, remainingPoints]  = distanceBasedClustering(X, maxDist, clusterRemainFlag)
    % X is the input data matrix with rows as points [x, y]
    % maxDist is the maximum allowable distance from the cluster center
    % if clusterRemainFlag is true than cluter the remainingPoint in
    % seperate clusters
    
    % input 
    if nargin < 3
        clusterRemainFlag = true;
    end
    
    % One Point Clutteriring
    if size(X,1) == 1
        clusters{1} = X;
        C = X;
        remainingPoints = [];
        return
    end
    
    % Estimate the number of clusters
    numClusters = estimateNumClusters(X, maxDist);
    numClusters = min(numClusters, size(X,1));
    % Initial clustering using k-means
    [idx, C] = kmeans(X, numClusters);
    
    % Create a cell array to hold the points in each cluster
    clusters = cell(0);
    usedPoints = false(size(X, 1), 1);
    
    % Assign points to clusters based on the maxDist constraint
    for i = 1:numClusters
        % Find points in this cluster
        clusterPoints = X(idx == i, :);
        
        % Calculate distances to the cluster center
        distances = sqrt(sum((clusterPoints - C(i, :)) .^ 2, 2));
        
        % Retain only those points within the specified distance
        validPoints = clusterPoints(distances <= maxDist, :);
        clusters{end+1} = validPoints; %#ok<AGROW>
        
        % Mark used points
        usedPoints(idx == i) = distances <= maxDist;
    end
    
    
    remainingPoints = X(~usedPoints, :);
    % Assign remaining points to their own clusters
    if clusterRemainFlag
        for i = 1:size(remainingPoints, 1)
            clusters{end+1} = remainingPoints(i, :); %#ok<AGROW>
            C(end+1,:) = remainingPoints(i, :);
        end
    end

%     figure;
% %     hold on;
%     colors = lines(length(clusters));
%     for i = 1:length(clusters)
%         scatter(clusters{i}(:,1), clusters{i}(:,2), 48, 'filled', 'MarkerFaceColor', colors(i,:));
%     end
% title('Custom Clustering with Maximum Distance Constraint');
% xlabel('X');
end % function


%%%%%%%
%% Estimate number of clusters
function numClusters = estimateNumClusters(X, maxDist)
    % Estimate the number of clusters needed based on the density and maxDist
    % This is a heuristic and can be adjusted based on your specific data
    numPoints = size(X, 1);
    area = pi * maxDist^2;
    density = numPoints / (range(X(:,1)) * range(X(:,2)));
    numClusters = ceil(numPoints / (density * area));
end

