function [] = PositionInEnvironmentCoordination(FindsStruct)
% THis Function clculates the findings of the bat in the simulation
% environmrnt coordinatiob (X-Y, not metric)
% Inputs:
%       FindsStruct: A Struct of the finding needed to be found
%           FindsStruct.Distance : the Distances from the Bat
%           FindsStruct.Alpha : THe angles from the bat
%           FindsStruct.IsAny : A flag of findings somthing
    
    RangeVec = FindsStruct.Distance;
    AngleVec = FindsStruct.Alpha;
    IsObsVec(nTime) = IsObs;
    xFindCurrentVec = xBati(nTime)+RangeVec.*cos(Teta(nTime)+AngleVec);
    yFindCurrentVec = yBati(nTime)+RangeVec.*sin(Teta(nTime)+AngleVec);
    tObsFind = [tObsFind, nTime];
    xObsFind = [xObsFind, xFindCurrentVec]; % update the find obsticle vectors
    yObsFind = [yObsFind, yFindCurrentVec];
end % function