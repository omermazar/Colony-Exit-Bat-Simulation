function [clusters, centers] = recursiveDistanceBasedClustering(X, maxDist)
    % Initialize clusters
    clusters = {};
    centers  = [];
    
    % Perform recursive clustering
    remainingPoints = X;
    numRemaining = size(remainingPoints,1);
    maxIter = 5;
    prevNum = 0;
    count = 0; % counter
    while (numRemaining ~= 0  || numRemaining > prevNum) && (count <= maxIter)
        prevNum = numRemaining;
        try 
            [newClusters, newCenters, remainingPoints] = distanceBasedClustering(remainingPoints, maxDist, false);
        catch
            warning(strcat('recurssive Clustering Error, count = ', num2str(count)))
        end

        % add current iteration to clusters
        clusters = [clusters,  newClusters];
        centers = [centers; newCenters];
        numRemaining = size(remainingPoints,1);
        count = count +1;
    end 
    % Add remaining points as individual clusters
    for i = 1:size(remainingPoints, 1)
        clusters{end+1} = remainingPoints(i, :); %#ok<AGROW>
    end

% figure;
% hold on;
% colors = lines(length(clusters));
% for i = 1:length(clusters)
%     scatter(clusters{i}(:,1), clusters{i}(:,2), 48, 'filled', 'MarkerFaceColor', colors(i,:));
%     plot(centers(i,1), centers(i,2), '*', 'MarkerSize',12, 'LineWidth', 2, 'Color', colors(i,:))
% 
% end
% title('Recursive Distance-Based Clustering');
% xlabel('X');
% ylabel('Y');

end % function





% %%%%%%%%%%%%%%%% 
% function [clusters, remainingPoints] = clusterPoints(points, maxDist, clusters)
%     if isempty(points)
%         remainingPoints = [];
%         return;
%     end
% 
%     % Estimate the number of clusters
%     numClusters = estimateNumClusters(points, maxDist);
% 
%     % Initial clustering using k-means
%     [idx, C] = kmeans(points, numClusters, 'MaxIter', 100);
% 
%     remainingPoints = [];
%     for i = 1:numClusters
%         % Find points in this cluster
%         clusterPoints = points(idx == i, :);
% 
%         % Calculate distances to the cluster center
%         distances = sqrt(sum((clusterPoints - C(i, :)) .^ 2, 2));
% 
%         % Retain only those points within the specified distance
%         validPoints = clusterPoints(distances <= maxDist, :);
%         clusters{end+1} = validPoints; %#ok<AGROW>
% 
%         % Collect remaining points for further clustering
%         remainingPoints = [remainingPoints; clusterPoints(distances > maxDist, :)]; %#ok<AGROW>
%     end
% 
%     % Recursively cluster the remaining points
%     if ~isempty(remainingPoints)
%         [clusters, remainingPoints] = clusterPoints(remainingPoints, maxDist, clusters);
%     end
% end

function numClusters = estimateNumClusters(X, maxDist)
    % Estimate the number of clusters needed based on the density and maxDist
    % This is a heuristic and can be adjusted based on your specific data
    numPoints = size(X, 1);
    area = pi * maxDist^2;
    density = numPoints / (range(X(:,1)) * range(X(:,2)));
    numClusters = ceil(numPoints / (density * area));
end

% Example usage
% Generate example data
