function [ixPoints] = filterFarSide(pointsInBeam, minDist, sortFlg)

% sort angles from right to left
[tetasOut, iSort] = sort(pointsInBeam.tetas);
distancesOut = pointsInBeam.distances(iSort);

% find the closest point
[~, imin] = min(distancesOut);
  
% Start from th mclosest point

% Filter out right-side
ixRightIn = iSort(imin :-1:1);
ixRight = funcFilter(pointsInBeam, ixRightIn, minDist);

% % Filter out left-side
ixLeftIn = iSort(imin : 1 : numel(tetasOut));
ixLeft = funcFilter(pointsInBeam, ixLeftIn, minDist);

ixPoints = [ixRight, ixLeft];
if nargin > 3 && sortFlg
    ixPoints= sort(ixPoints);
end

end % main
%% 
function [ixClose] = funcFilter(pointsInBeam, ixSide, minDist)

% Init
xin = pointsInBeam.xin;
yin = pointsInBeam.yin;
distances = pointsInBeam.distances;
ixClose = ixSide;
k=1;

% start the loop
while k < numel(ixClose)
    ix = ixClose(k);
    ixCheck = ixClose(k+1:end);
    ixN = find(sqrt((xin(ixCheck)-xin(ix)).^2 + (yin(ixCheck)-yin(ix)).^2) <= minDist, 1);
    %     plot(xin(ix), yin(ix), 'k*')
    %     plot(xin(ixCheck), yin(ixCheck), 'ok')

    if ixN > 1
        % find the next minimal distance
        ixJump = find(distances(ixClose((k+1) : (k+ixN-1))) < distances(ix) , 1);
        if isempty(ixJump) % ther are no closer point between remove the inbetween point
            ixClose((k+1):(k+ixN-1)) = [];
        else % if
            ixClose((k):(k+ixJump-1)) = [];
            k=k-1; % check again the Jump Point
        end
    end % if ixN
    % go to next point
    k = k+1;
    if k > numel(xin)
        warning('while loop exceed vector diemnsions')
        break
    end
end % while

%%% test the Last point, and remove it if needed
if numel(ixClose) > 1
    minEdgeResolution = 50*1.02;
    distLast = sqrt((xin(ixClose(end))-xin(ixClose(end-1))).^2 + (yin(ixClose(end))-yin(ixClose(end-1))).^2);
    if distLast >= minDist && distLast <= minEdgeResolution
        ixClose(end) = [];
    end % if distLast >= minDist && distLast <= minEdgeResolution
end %if ~isempty(ixClose)
end % fuction

