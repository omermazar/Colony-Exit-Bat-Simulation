function [PositionOfOtherBats] = FindPosFromBatsStrct(BAT, MyBatNum, nTime, NumberOfBats)

% this funcion return the posions of the other bats in the environment
% output:
%   PositionOfOtherBats ((NumOfBats-1)x5 matrix) - each raw with fields: 
%           [BatNumber, xBat(nTime), yBat(nTime), BatTeta(nTime), Batvelocity(nTime)]
%
% inputs:
%  BAT - the strcuture of the BATS
%  MyBatNum - the current BAT - not intersting
%  Ntime- the Time intersted
 
if NumberOfBats >1
    PositionOfOtherBats = zeros(NumberOfBats-1,5);
    RawNum=1;
    for kk= 1:NumberOfBats
         if kk ~= MyBatNum
            ObatNumber =  kk;
            OBatxPos = BAT(kk).xBati(nTime);
            OBatyPos = BAT(kk).yBati(nTime);
            OBatTeta = BAT(kk).Teta(nTime);
            OBatVel = BAT(kk).BatVelocity(nTime);
            PositionOfOtherBats(RawNum,:)  = ...
                [ ObatNumber, OBatxPos, OBatyPos, OBatTeta, OBatVel ];
            RawNum = RawNum +1;
         end % if kk ~= MyBatNum
    end % for kk= 1:NumberOfBats
else % if NumberOfBats >1
    PositionOfOtherBats = [];
end % if NumberOfBats >1


