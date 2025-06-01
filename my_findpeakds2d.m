function [pks, row, column] =my_findpeakds2d(data, varargin) 

data_tr = data';
[pks1, locs1, w1, p1] = findpeaks(data(:)) ; % peaks along x
[~, locs2, ~, ~] = findpeaks(data_tr(:)); % peaks along y

data_size = size(data); % Gets matrix dimensions
[row2, col2] = ind2sub(size(data_tr), locs2); % Converts back to 2D indices
locs2 = sub2ind(data_size, col2, row2); % Swaps rows and columns and translates back to 1D indices

[ipks, i1, i2] = intersect(locs1, locs2); % Finds common peak position
[row, column] = ind2sub(data_size, ipks); % to 2D indices
pks = pks1(i1);


