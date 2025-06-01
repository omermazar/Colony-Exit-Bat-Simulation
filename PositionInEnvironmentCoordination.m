function [xFinds, yFinds] = PositionInEnvironmentCoordination(xBat, yBat, ViewAngle, FindsStruct, varargin)
% THis Function clculates the findings of the bat in the simulation
% environmrnt coordinatiob (X-Y, not metric)
% Inputs:
%       xBat, yBat- the current Position of the bat
%       ViewAngle - the angle the bat is looking (flying) to 
%       FindsStruct: A Struct of the findings needed to be found (obs or preys), can
%       include vector of all the findings
%           FindsStruct.Distance : the Distances from the Bat
%           FindsStruct.Alpha : THe angles from the bat
%           FindsStruct.IsAny : A flag of findings somthing
%      If THe struct is for PREY- use: 
%           PositionInEnvironmentCoordination(xBati(nTime),yBati(nTime),Teta(nTime), PreyStruct, ...
%                        NumOfPreyDetected, DetectedPreys)
% OutPuts:
%       xFinds, yFinds = vectors of the coordination found 

NumberOfInputs= nargin;
% % % if NumberOfInputs == 5
% % %     NumOfEchos =  varargin{1};
% % %     DetectedPreys = varargin{2};
% % %     FindsRangeVec = FindsStruct.Distances;
% % %     FindsAngleVec = FindsStruct.TargetAngle;
% % % else
% % %     FindsRangeVec = FindsStruct.Distances;
% % %     FindsAngleVec = FindsStruct.TargetAngle;
% % % end %if NumberOfInputs == 6
 
FindsRangeVec = FindsStruct.Distances;
FindsAngleVec = FindsStruct.TargetAngle;

xFinds = xBat + FindsRangeVec.*cos(ViewAngle + FindsAngleVec);
yFinds = yBat + FindsRangeVec.*sin(ViewAngle + FindsAngleVec);

end % function