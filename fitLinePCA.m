function [a, b] = fitLinePCA(X, Y)
    % Combine X and Y into a single data matrix
    data = [X(:) Y(:)];
    
    % Center the data
    data_mean = mean(data, 1);
    data_centered = data - data_mean;
    
    % Perform PCA
    [coeff, ~, ~] = pca(data_centered);
    
    % The direction of the line is given by the first principal component
    direction = coeff(:, 1);
    
    % Calculate the slope (a) and intercept (b) of the line
    a = direction(2) / direction(1);
    b = data_mean(2) - a * data_mean(1);
end