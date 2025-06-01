function pairs = findClosePairs(matrix, minDistThreshold)
    % %%%% New June2024
    % Find pairs of rows in a matrix with Euclidean distance less than minDistThreshold
    %
    % Inputs:
    %   matrix - a 2D matrix where each row is a point
    %   minDistThreshold - the distance threshold
    %
    % Outputs:
    %   pairs - an Nx2 matrix of row indices that are closer than minDistThreshold

    % Calculate pairwise Euclidean distances
    distances = pdist2(matrix, matrix);

    % Find pairs of rows with distances below the threshold
    [rowIdx, colIdx] = find(tril(distances,-1) < minDistThreshold)
    [rowIdx, colIdx] = find( distances < 0.1);
    
    ixDiag = rowIdx ~= colIdx;
    % Combine indices into pairs
    pairs = [rowIdx(ixDiag), colIdx(ixDiag)];
end



matrix = [[BAT.ObsFindsStruct(minPulse:AnalyzedPulseNum-1).xFinds]', [BAT.ObsFindsStruct(minPulse:AnalyzedPulseNum-1).yFinds]'];
% Find unique rows and their indices
[uniqueRows, ~, uniqueIndices] = unique(matrix, 'rows');

% Count occurrences of each unique row
occurrences = histcounts(uniqueIndices, unique(uniqueIndices));
 numUnique = max(uniqueIndices)
occurrences2 = histcounts(uniqueIndices, 0.5:1:(numUnique + 0.5)); % histcounts(uniqueIndices, 'BinMethod', 'integers')
% Find rows that occur more than once
duplicateRowIndices = find(occurrences2 > 1);
duplicateRows = uniqueRows(duplicateRowIndices, :);

% Display the duplicate rows
disp('Rows that occur more than once:');
disp(duplicateRows);