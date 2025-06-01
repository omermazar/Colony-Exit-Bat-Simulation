function [out_timesA, out_pksA, detected_ix, time_diff] =  ...
    detections_compare(in_timesRef, in_timesA, in_pksA, time_tol)
% Compare between two sets of detection, by time 
%%% NOV2021

%init
detected_ix = false(size(in_timesRef));
% protect against emty struct
% if isempty(in_timesRef)
%     in_timesRef = inf;
% end

% init
out_timesA = [];
time_diff  = [];
out_pksA   = [];

% check if none is detected
if ~isempty(in_timesA) & ~isempty(in_timesRef) 
    [k,dist] = dsearchn(in_timesRef', in_timesA');

    locB = find(dist < time_tol);
    
    %% find unique refs only)
    k_found = k(locB);
    detected_ix(k_found) = true;
    if ~isempty(k_found) % numel(k_found) ~= numel(unique(k_found))
        dist_found = dist(locB);
        g          = findgroups(k_found);
        %         [~,~,g] = unique(k_found,'stable');
        % number of points in each group
        pointsInGroup = splitapply(@numel, k_found, g);
        % find the minmum distance to each reference point
        [dist_min, minIx] = splitapply(@min, dist_found, g);
        % find the index in the origina of the mimal distances
        nn = cumsum([0;pointsInGroup(1:end-1)])+minIx;
        % % %        locB = ismember(dist, dist_min);
        locA = locB(nn);
        %%% output
        out_timesA = in_timesA(locA); % in_timesA(locB);
        time_diff  = out_timesA - in_timesRef(k(locA));% out_timesA - in_timesRef(k(locB));
        if nargout > 1
            out_pksA = in_pksA(locA);
        end
   
    end % if isempty(k_found) % if numel(k_found) ~= numel(unique(k_found))
    

else % if ~isempy(in_timesA)
     %% output
    out_timesA = [];
    time_diff  = [];
    out_pksA   = [];
end % if ~isempy(in_timesA)
