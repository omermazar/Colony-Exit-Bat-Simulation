function [ManueverCmdStruct] = keepDistancefromObs(ManueverCmdStruct, BAT, AnalyzedPulseNum,  ...
    sortedAngles, sortedDistancesByAngles, sortedXEst, sortedYEst, ...
    distToKeepFromWall, MinDistanceAllowed, DistanceToBuzz)

%%% New June2024 - not used
% parameters
inFrontAngleLimit = pi/4; % 45 dergrees

%% the postion of the bat
nTime = BAT.TransmittedPulsesStruct(AnalyzedPulseNum).StartPulseTime;
xBat  = BAT.xBati(nTime);
yBat  = BAT.yBati(nTime);
tetaB = BAT.Teta(nTime);

%% the refence planned point to navgate to
refDist = ManueverCmdStruct.Dist2HuntedPrey;
requiredDirection = ManueverCmdStruct.Angle2HuntedPrey ;
plannedPoint = [xBat + refDist*cos(tetaB+requiredDirection), yBat + refDist*sin(tetaB+requiredDirection)];

%% Look for the closest obstable infront of the bat and adjust flight path
ixInfornt = abs(sortedAngles) <= inFrontAngleLimit;
% the closest distance to the planned flight path
x0 = sortedXEst(ixInfornt);
y0 = sortedYEst(ixInfornt);
x1 = xBat; y1= yBat;
x2 = plannedPoint(1); y2 = plannedPoint(2);
% Calculate the distance using the formula
numerator = abs((y2-y1).*x0 - (x2-x1).*y0 + x2.*y1 - y2.*x1);
denominator = sqrt((y2-y1).^2 + (x2-x1).^2);
d = numerator ./ denominator;
minDist   = min(d);

% if too close slow down and correct the angle
if minDist < distToKeepFromWall
    ManueverCmdStruct.ManAccelCommand = 'SlowDown';    
    requiredDirection = requiredDirection + sign(requiredDirection)*pi/6;
    ManueverCmdStruct.Angle2HuntedPrey = requiredDirection;
end

% Buzz
if minDist < DistanceToBuzz
end


