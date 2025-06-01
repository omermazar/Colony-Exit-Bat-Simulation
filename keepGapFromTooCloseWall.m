%% %%%%% Keep Gap from Wall
function requiredDirection = keepGapFromTooCloseWall(requiredDirection, checkAngles, checkDistances, distToKeepFromWall)

DistfromWall = abs(sin(requiredDirection - checkAngles) ) .* checkDistances;
[minD, minIxx] = min(DistfromWall);
angleFromWall = 10/180*pi;
if minD < distToKeepFromWall/3
    requiredDirection = requiredDirection + sign(requiredDirection- checkAngles(minIxx))*angleFromWall;
end

end % function