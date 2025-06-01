function [requiredDirection, avoidCmd, minDist, minAngle, objInd] = avoidObstacleInfornt( ManueverType, requiredDirection, redquiredDist, ...
                    xBat, yBat, tetaB, sortedAngles, sortedDistancesByAngles, sortedXEst, sortedYEst, ...
                    inFrontAngleLimit, distToKeepFromWall, MinDistanceAllowed, MinManueverRadius)

% this function gets the planned manuver of the bats and the detected
% obstacles  and returns the reuired diriecton to avoid Obstacles
% avoidCmd = "none" | "avoid" | "crushAlert"
%  MinManueverRadius = CurrentVelocity.^2/MaxAccel

%% Input
if isempty(redquiredDist)
    redquiredDist = 1.5 / 0.01;
end

avoidCmd = "none";
correcrtionAngle = pi/15; 
distDiffToDecide = distToKeepFromWall; %0.8 / 0.01;

%% Look for the closest obstable infront of the bat and adjust flight path
ixInfront = abs(sortedAngles) <= inFrontAngleLimit;
% the closest distance to the planned flight path
xObj    = sortedXEst(ixInfront);
yObj    = sortedYEst(ixInfront);
angObj  = sortedAngles(ixInfront);
distObj = sortedDistancesByAngles(ixInfront);

%%% caulcate distances between all point to the planned path
% the palnned path
x1 = xBat; y1= yBat;
% x2 = plannedPoint(1); y2 = plannedPoint(2);
x2 = x1 + redquiredDist*cos(requiredDirection+tetaB);
y2 = y1 + redquiredDist*sin(requiredDirection+tetaB);
% xx = linspace(x1,x2,10);
% yy = linspace(y1,y2,10);
% the distances
% allDists = min(pdist2([xObj',yObj'],[xx',yy']),[],2);
% [minDist, mi] = min(allDists);
% Calculate the distance using the formula

numerator = abs((y2-y1).*xObj - (x2-x1).*yObj+ x2.*y1 - y2.*x1);
denominator = sqrt((y2-y1).^2 + (x2-x1).^2);
d = numerator ./ denominator;
[~, mi]   = min(d);
minAngle = angObj(mi);
objInd = mi;

% find closest obstacle in each side
indRight = find(angObj < requiredDirection); % negative, 4th qurter) 
indLeft  = find(angObj >= requiredDirection); % positive, 1st quarter
if any(indRight)
    [minDistRight, miRight] = min(d(indRight));
else
    minDistRight = inf;
end
if any(indLeft)
    [minDistLeft,  miLeft] = min(d(indLeft));
else
    minDistLeft = inf;
end
[minDist, miRL] = min([minDistRight, minDistLeft]); % miRL = 1 if right obstacle is closer to the line, 2 if left

fixAngle = 0;
   
%%%5 check if the point is above or below the line - OLD
% if minAngle >= requiredDirection % abobe (positive)
%     fixAngle = -correcrtionAngle;
% else % below
%     fixAngle = correcrtionAngle;
% end

%%% the reaction is based on manueverType
switch ManueverType
    % if too close slow down and correct the angle
    case 'FollowWall'
        distToKeep = 0.75*distToKeepFromWall;

    case 'ExitMan'
        distToKeep = 0.5*distToKeepFromWall;
end % switch

%%%% check for crush in front
if minDist < MinDistanceAllowed 
    avoidCmd = "crushAlert";
    if minDistLeft < MinDistanceAllowed && minDistRight < MinDistanceAllowed
         % too close obstacles in both sides
        % avoid the obstacle that is closer to the current point
        if distObj(indLeft(miLeft)) < (distObj(indRight(miRight)) - distDiffToDecide)
            % turn right
            fixAngle = -correcrtionAngle;
        elseif distObj(indRight(miRight)) < (distObj(indLeft(miLeft)) - distDiffToDecide)
            % turn Left (easy)
            fixAngle = correcrtionAngle;
        else %  Turn to the 'open side" 
            if max(abs(angObj(indLeft))) <= max(abs(angObj(indRight)))
               % turn Left (easy)
                fixAngle = correcrtionAngle; 
            else
                % turn right
                fixAngle = -correcrtionAngle;
            end
            % stay ahed
            % fixAngle = 0;
        end
    elseif miRL == 1 % Right is the dangere, turn Left 
        fixAngle = 2*correcrtionAngle;
    else % Turn right
        fixAngle = -2*correcrtionAngle;
    end

%%%% Check for keepDist
elseif minDist < distToKeep 
    avoidCmd = "avoid";
    if minDistLeft < distToKeep && minDistRight < distToKeep
        % obstacles in both sides - stay ahed
        fixAngle = 0;
    elseif miRL == 1 % Right is the dangere, turn Left 
        fixAngle = correcrtionAngle;
    else % Turn right
        fixAngle = -correcrtionAngle;
    end
%%%%no response is needed
else
       %no response
end % if minDist < 2*MinDistanceAllowed 

requiredDirection = requiredDirection + fixAngle;

%%%% Check if the required manuver is possible in the cuurent celocity
if requiredDirection <= 0 % if turn right
    ixConsider = angObj < 0  & abs(angObj) <= abs(requiredDirection);
else
    ixConsider = angObj > 0  & angObj <= abs(requiredDirection);
end
% check all obsacle beteen cuurent directoin and the reuired direction
if min(distObj(ixConsider)) <= 0.5*MinManueverRadius
     if ~strcmp(avoidCmd,'crushAlert')
         avoidCmd = "avoid";
     end
elseif min(distObj(ixConsider)) <= MinManueverRadius
    avoidCmd = "crushAlert";
end





   % calcualte the line
    % m = (y2 - y1) / (x2 - x1);
    % b = y1 - m * x1;
    % % chek
    % yCheck = m*sortedXEst(mi) + b;
    % if sortedYEst(mi) >= yCheck % abobe (positive)
    %     fixAngle = -correcrtionAngle;
    % else % below
    %     fixAngle = correcrtionAngle;
    % end
