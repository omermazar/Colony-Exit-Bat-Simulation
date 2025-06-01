function [TurnAngleDesired, DirectionCommand, AccelCommand, ManueverPower]= ObsDirectionDecision(FoundObstacles,AllParams,BeamWidth,CurrentVelocity,PrevDirectionCmd)
% This function finds the desired bat flight's direction with regarded to the Obstacles found in the arena
% Thisfunction is called only if IsObs==1 or ~isempty(FoundObsticles)
% inputs:
% FoundObstacles:structure of the obstacle found (distance, angle from ...
% ... the bat, Isobs flag)
% BeamWidth: The detection Sonar Beam-Width
% Simparameters - defult params
% CurrentVelocity - the Bat's Cuurent Velocity

% Outputs:
% TurnAngleDesired - the angle to turn
% DirectionCommand - 'Left' / 'Right'
% AccelCommand     - 'SlowDown' / 'None' / 'Accelerate'
% ManueverPower    - 'RegularManuever' / 'CrushAvoidance'
[TargetAngle, IndexSorting] = sort(FoundObstacles.TargetAngle); 
FoundDistances  = FoundObstacles.Distances(IndexSorting);

    
delta_angle     = zeros(1,length(TargetAngle)+1); %initiate
AccelCommand    = 'None'; %The default command for Acceleration. in case of 'slowdown' it will be specified in the code
ManueverPower   = 'RegularManuever';%The default command for ManeuverPower. in case of 'CrushAvoidance' it will be specified in the code
for i=1:1:length(TargetAngle)+1
    if i==1
        delta_angle(i) = abs(-BeamWidth/2-TargetAngle(i)); %dealing with the first element
    elseif i==length(TargetAngle)+1
        delta_angle(i) = abs(BeamWidth/2-TargetAngle(end)); %dealing with the last element
    else
        delta_angle(i) = TargetAngle(i)-TargetAngle(i-1);
    end
end

NoObsZone = find(abs(delta_angle)>=10*pi/180); %find where there are No obstacles within 10 degrees between them
ThetaZone = zeros(length(NoObsZone),2); %initiate
if  isempty(NoObsZone) %all the FOV is coverd with Obstacles
    TurnAngleDesired = pi/2;
    [IsRightManeuverPossibleFlag, RightDist] = CheckDirectionManeuver('Right', FoundDistances, TargetAngle, CurrentVelocity, AllParams,-TurnAngleDesired);
    [IsLeftManeuverPossibleFlag, LeftDist]   = CheckDirectionManeuver('Left',  FoundDistances, TargetAngle, CurrentVelocity, AllParams,TurnAngleDesired);
    % Decide which direction
    if IsRightManeuverPossibleFlag %if we can turn right 
        if ~IsLeftManeuverPossibleFlag % but not left - choose right
            DirectionCommand = 'Right';
        else  % if we can turn left and right - check distance and turn to the faraway side
            if RightDist<100 && LeftDist<100 % if both distances are very close to the wall, continue with previous command
                DirectionCommand = 'Right';%'PrevDirectionCmd';
            elseif max(RightDist)<= max(LeftDist)
                DirectionCommand = 'Left';
            else 
                DirectionCommand = 'Right';
            end
        end
    else %if we can't turn right
        if IsLeftManeuverPossibleFlag %if we can turn left
            DirectionCommand    = 'Left';
        else % otherwise (we can't turn left or right) - slow downd and turn to the faraway side
            AccelCommand        = 'SlowDown';
            ManueverPower       = 'CrushAvoidance';
%             DirectionCommand    = 'Right';
            if RightDist<100 && LeftDist<100 % if both distances are very close to the wall, continue with previous command
                DirectionCommand =  'Right'; %'PrevDirectionCmd';
            elseif max(RightDist)<= max(LeftDist)
                DirectionCommand = 'Left';
            else 
                DirectionCommand = 'Right';
            end
        end 
    end 
else %we have fields without obstacles
    for i=1:1:length(NoObsZone)
        if NoObsZone(i)==1 % at the left edge of the FOV
            ThetaZone(i,:) = [-BeamWidth/2,TargetAngle(NoObsZone(1))];
        elseif NoObsZone(i)== length(delta_angle) %at the right edge of the FOV
            ThetaZone (i,:)= [TargetAngle(NoObsZone(end)-1),BeamWidth/2];
        else %in the middle range of the FOV
            ind            = NoObsZone(i);
            ThetaZone(i,:) = [TargetAngle(ind-1),TargetAngle(ind)];
        end
    end
    MidThetaZone         = (ThetaZone(:,1)+ ThetaZone(:,2))./2;
    Range                = zeros(length(MidThetaZone),1); %initiate
%     newVec_MidTheta      = zeros(length(MidThetaZone),1); %initiate
    [Max,Indx]           = max(delta_angle(NoObsZone));
    theta_desired        = MidThetaZone(Indx);
    if ThetaZone(Indx,1)<=-7.5*pi/180 && ThetaZone(Indx,2)>=7.5*pi/180 %if in the flight path direction ±5 deg there are no obstacles - continue straight
            TurnAngleDesired = 2*pi/180*rand(1)-2*pi/180*rand(1);  %random degree between -2 to 2
            DirectionCommand = 'None';
    else % Check if it possible to turn
        TurnAngleDesired    = theta_desired;%+2*pi/180*(-1+2*rand(1));
        if TurnAngleDesired < 0
           DirectionCommand = 'Left';
        else 
           DirectionCommand = 'Right';
        end
        [IsDirectionPossibleFlag, Range(1)]    = CheckDirectionManeuver(DirectionCommand, FoundDistances, TargetAngle, CurrentVelocity, AllParams,TurnAngleDesired);
        if ~IsDirectionPossibleFlag %if we cannot turn to the max Zone
            for i =2:1:length(MidThetaZone)
                a                   = find(delta_angle(NoObsZone)<Max);
                [Max,Indx]          = max(delta_angle(a)); 
                theta_desired       = MidThetaZone(Indx);
                if ThetaZone(Indx,1)<=-7.5*pi/180 && ThetaZone(Indx,2)>=7.5*pi/180 %if in the flight path direction ±5 deg there are no obstacles - continue straight
                    TurnAngleDesired = 2*pi/180*rand(1)-2*pi/180*rand(1);  %random degree between -2 to 2
                    DirectionCommand = 'None';
                    break
                else % Check if it possible to turn
                    if theta_desired < 0
                       DirectionCommand = 'Left';
%                        sign_for_angle   = -1;
                    else 
                       DirectionCommand = 'Right';
%                        sign_for_angle   = 1;
                    end
                    TurnAngleDesired                    = theta_desired;%+2*pi/180*(-1+2*rand(1));
                    [IsDirectionPossibleFlag, Range(i)] = CheckDirectionManeuver(DirectionCommand, FoundDistances, TargetAngle, CurrentVelocity, AllParams,TurnAngleDesired);                  
    %                 [~,ind3]                          = min(abs(ThetaZone(Indx,:)));
    %                 TurnAngleDesired                  = ThetaZone(Indx,ind3)+sign_for_angle*10*pi/180;%+2*pi/180*(-1+2*rand(1));
                    if ~IsDirectionPossibleFlag         %the turn is not possible
    %                    newVec_MidTheta(i)             = MidThetaZone(Indx);
                       MidThetaZone(Indx)               = 1000;
                    else
                        break
                    end
                    if i==length(MidThetaZone)-1&& MidThetaZone(Indx) == 1000 %if we didnt find any possible zone to turn without Crushing the walls
                        [~,Index2]          = max(Range);
                        TurnAngleDesired    = MidThetaZone(Index2);%+2*pi/180*(-1+2*rand(1));
                        AccelCommand        = 'SlowDown';
                        ManueverPower       = 'CrushAvoidance';          
                        if TurnAngleDesired < 0  
                           DirectionCommand = 'Left';
                        else
                           DirectionCommand = 'Right';
                        end
                    end
                end
            end            
        end
    end
end
