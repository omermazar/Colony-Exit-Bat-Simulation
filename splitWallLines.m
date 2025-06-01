function [newEdges, newlines, newPolys] = splitWallLines(sortedXest, sortedYest, min_points_per_line, plotFlag)
% refLine = [m,b]; %y = m*x + b

% input 
x = sortedXest;
y = sortedYest;

% Sample data: replace with your (x, y) data
% x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
% y = [1, 2.1, 3.2, 6, 7, 8.1, 13, 13.9, 15, 15.5]; % Example with noise

% Parameters
window_size = 5; % Size of the window for linear fitting
slide_size = 5; % Slide by 4 pointsslope_threshold = 1; % Slope change threshold
slope_threshold = 1; % Slope change threshold
angle_threshold = pi/4; % 45 degrees between walls
intercept_threshold = 2; 100; % Intercept change threshold
% min_points_per_line = 3; % Minimum points in a line segment

% Initialize variables
numPoints = length(x);
newlines = {};
indices = {}; % Store indices for each segment
newPolys = [];

current_line = [x(1:window_size)', y(1:window_size)'];
current_indices = 1:window_size;

% Fit a line to each window
i = 1;
while i <= numPoints % i <= (numPoints - window_size + 1)
    % Current window
    endInd = min([numPoints, i + window_size - 1]); % the last window may be shorter
    x_window = x(i:endInd); % x(i:(i + window_size - 1))
    y_window = y(i:endInd); % y(i:(i + window_size - 1))
    
    % if the last segemnet is less than minLength - add it to the prevoius segement
    % Add new points from the current window to the current line
    % Avoid adding duplicate points
    if numel(x_window) < min_points_per_line
        % new_points_start = current_line(end, 1) + 1; % Start after the last point in current_line
        % new_points_index = find(x_window >= new_points_start);
        % current_line = [current_line; x_window(new_points_index)', y_window(new_points_index)'];
        % current_indices = [current_indices, i + new_points_index - 1]; % Update indices
        
        new_indices = setdiff(i:endInd, current_indices);
        % trim current_indices to the number of point in lines
        new_indices = new_indices(new_indices<=endInd);
        current_line = [current_line; x(new_indices)', y(new_indices)'];
        current_indices = [current_indices, new_indices]; % Update indices

        % exit the while loop
        i = inf;
        continue
    end

    % Linear regression for the current window
    p = fit_near_vertical_line(x_window', y_window'); % polyfit(x_window, y_window, 1)
    slope_current = p(1);
    intercept_current = p(2);

    % Linear regression for the previous window (if not the first iteration)
    if i > 1
        x_prev_window = x((i - slide_size):(i - slide_size + window_size - 1));
        y_prev_window = y((i - slide_size):(i - slide_size + window_size - 1));
        p_prev = fit_near_vertical_line(x_prev_window', y_prev_window'); % polyfit(x_prev_window, y_prev_window, 1)
        slope_prev = p_prev(1);
        intercept_prev = p_prev(2);

        % Check if the difference exceeds the threshold
        % Calculate angle between lines using arctan
        angle_change = abs(atan((slope_current - slope_prev) / (1 + slope_current * slope_prev)));
        
        % Check if the difference exceeds the angle and intercept thresholds
        if angle_change > angle_threshold % || abs(intercept_current - intercept_prev) > intercept_threshold
            % If current line has at least min_points_per_line, save it
            if size(current_line, 1) >= min_points_per_line
                % calculaete the linear regression for the full line
                pLine = fit_near_vertical_line(current_line(:,1), current_line(:,2)); % polyfit(current_line(:,1), current_line(:,2), 1)

                % save the current line
                newlines{end+1} = current_line; % Save the current line
                indices{end+1} = current_indices;
                % save he regerssion
                newPolys = [newPolys; pLine];
            end
            % Start a new line
            current_line = [x_window', y_window'];
            current_indices = i:(i + size(current_line,1) - 1);

        else
            % Add new points from the current window to the current line
            % Avoid adding duplicate points
            % new_points_start = current_line(end, 1) + 1; % Start after the last point in current_line
            % new_points_index = find(x_window >= new_points_start);
            % current_line = [current_line; x_window(new_points_index)', y_window(new_points_index)'];
            % current_indices = [current_indices, i + new_points_index - 1]; % Update indices

            % Add all points from the current window to the current line
            % Append new points and update indices without assuming order
            new_indices = setdiff(i:endInd, current_indices);
            current_line = [current_line; x(new_indices)', y(new_indices)'];
            current_indices = [current_indices, new_indices]; % Update indices
        end
    end
    
    % Move the window by the slide size
    i = i + slide_size;
end % while

% Add the last line if it meets the condition
if size(current_line, 1) >= min_points_per_line
    newlines{end+1} = current_line;
    indices{end+1} = current_indices;
    % calculaete the linear regression for the full line
    pLine = fit_near_vertical_line(current_line(:,1), current_line(:,2)); % polyfit(current_line(:,1), current_line(:,2), 1)
    % save he regerssion
    newPolys = [newPolys; pLine];
end

% output
newEdges = nan(numel(newlines),2);
for k = 1:numel(newlines)
    newEdges(k,:) = [indices{k}(1), indices{k}(end)];
end

% Display results
if plotFlag
    fprintf('Number of lines found: %d\n', length(newlines));
    for k = 1:length(newlines)
        fprintf('Line %d:\n', k);
        disp(newlines{k});
    end

    % Plotting for visualization
    figure;
    hold on;
    plot(x, y, 'ko-', 'DisplayName', 'Original Data');
    colors = lines(length(newlines));

    for k = 1:length(newlines)
        line_segment = newlines{k};
        plot(line_segment(:,1), line_segment(:,2), '-', 'Color', colors(k,:), 'LineWidth', 2, 'DisplayName', ['Line ' num2str(k)]);
    end

    xlabel('X');
    ylabel('Y');
    title('Detected Line Segments');
    legend show;
    grid on;
    % hold off;
end