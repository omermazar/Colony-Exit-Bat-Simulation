function [DistVec, AngleVec, varargout] = FindDistAngle(xBat, yBat, BatDirection, xTarget, yTarget, varargin)
% this function caulculates the distance and angle between the vector points of the
% target(xTarget, yTarget) and the vector of the bat movement(xBat, yBat) and
% diercton of the movement (BatDirection)
    
    NofInputs = nargin;
    NofOutputs = nargout;
    TargetDirection = zeros(size(xTarget)); % default init
    if NofInputs == 6
        TargetDirection = varargin{1};
    end % if NofInputs == 6
   
    BatDirectionVec = [cos(BatDirection); sin(BatDirection)];
    DistVec = sqrt ( (yBat - yTarget).^2 + (xBat - xTarget).^2 );
    DeltaVec = [(xTarget- xBat) ; (yTarget - yBat)];
    AngleVec1 = acos( dot(BatDirectionVec,DeltaVec)./ DistVec); % might be plus or minus : aambiguity
        
    % Solving the Ambiguity of Cos - Finding the right angle x
%     ERR = 0.001*length(xBat);
    True_Angle= AngleVec1 + BatDirection;
    R1 = findHipDist(xBat,yBat,DistVec,True_Angle,xTarget,yTarget);
    AngleVec = zeros(1,length(R1)); % INIT
    for nn= 1: length(R1)
        ERR = 0.01*DistVec(nn); % the error is compared to the distance
        if R1 <= ERR
            AngleVec(nn) = AngleVec1(nn);
        else
            AngleVec(nn) = -AngleVec1(nn);
        end % if R1
    end % nn
    % Realative Direction
    if NofOutputs == 3
        B2TRelativeDirection = TargetDirection - BatDirection;
        
        % [-pi,pi]
%         B2TRelativeDirection = pi*sawtooth(B2TRelativeDirection-pi);
        B2TRelativeDirection = B2TRelativeDirection -(2*pi)*floor((B2TRelativeDirection+pi)/(2*pi));
        
        varargout{1} = B2TRelativeDirection ;
    end % if NofOutputs == 3

%% new - Cave Exit 20Nov22
%%%% deal with nans (when the bats exit the cave
ixNan           = isnan(DistVec );
DistVec(ixNan)  = inf;
AngleVec(ixNan) = pi;

end % function FindDistAngle

%%
%%%%
function [R] = findHipDist(xCurrent,yCurrent,PreyDist,PreyAbsAngle,xP,yP)
    XPreyHip =  xCurrent + PreyDist.*cos(PreyAbsAngle);
    YPreyHip =  yCurrent + PreyDist.*sin(PreyAbsAngle);
    R = sqrt((XPreyHip-xP).^2 + (YPreyHip-yP).^2 );
end % FindHipDist
