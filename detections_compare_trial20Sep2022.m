function [out_timesA, out_pksA, detected_ix, time_diff] =  ...
    detections_compare(in_timesRef, in_timesA, in_pksA, time_tol)
% Compare between two sets of detection, by time 
%%% NOV2021

%init
detected_ix = false(size(in_timesRef));
% protect against emty struct
if isempty(in_timesRef)
    in_timesRef = inf;
end

% check if none is detected
if ~isempty(in_timesA)
    [k,dist] = dsearchn(in_timesRef', in_timesA');

    locB = dist < time_tol;
    
    %% find unique refs only)
    k_found = k(locB);
    detected_ix(k_found) = true;
    if numel(k_found) ~= numel(unique(k_found))
        dist_found = dist(locB);
        uk = unique(k_found);
        N = numel(uk);
        ixB = nan(1,N);
        for ii = 1:N
            kRef = uk(ii);
            kIx = find(k_found == kRef);
            [~, minIx] = min(dist_found(kIx));
            ixB(ii) = kIx(minIx);    
        end % for kRef
%         g = findgroups(k_found);
%         dist_min = splitapply(@min, dist_found, g   );
%         locB = ismember(dist, dist_min);
    end % if numel(k_found) ~= numel(unique(k_found))
    
    %% output
%     out_timesA = in_timesA(locB);
%     time_diff  = out_timesA - in_timesRef(k(locB));

    out_timesA = in_timesA(ixB);
    time_diff  = out_timesA - in_timesRef(uk);

    if nargout > 1
        out_pksA = in_pksA(locB);
    end

else % if ~isempy(in_timesA)
     %% output
    out_timesA = [];
    time_diff  = [];
    out_pksA = [];
end % if ~isempy(in_timesA)
