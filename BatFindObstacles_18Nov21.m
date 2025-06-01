
function [ ObsStruct ] = ...
        BatFindObstacles(xCurrent,yCurrent,AngleView,BeamWidth, ...
        DetectionDistance, Terrain)
    
% This funtion will find if there are obsticles in the detection range...
% and return the range and angle from the current point
% Outputs:
    % ObsStruct - structure of the obstacle found:
    %       ObsStruct.Distances is the a vector of Distances between the bat and the
    %       obstacles found
    %       ObsStruct.TargetAngle are the angles between the bat flight direction and the obstacle
    %       ObsStruct.IsObs -1 if there was an obstacle Found    
% Inputs: 
    % xCurrent,yCurrent, - position and direction of the bat
    % AngleView - The Direction of the Bat Head (for now equals to his direction of
    % movement)
    % BeamWidth, DetectionDistance - The detection Sonar Beam-Width and
    % distance in the terrain coordinations system (no-units)
    % Terrain = the environment to check
    % PreysPositions = [xPrey,yPrey,PreyTeta, PreyVelictu] X Number Of Preys 
    
    %%% INIT %%%
    IsObs=0; 
    dt=0.1; % the step size for the angle in radians
  
    %%%% find Obstacles %%%
    
    LIMITS = size(Terrain);
    LIMITS_1 = LIMITS(1);
    LIMITS_2 = LIMITS(2);
    ObsDistance=zeros(1,ceil(DetectionDistance/dt)); % INIT
    ObsAlpha=zeros(1,ceil(BeamWidth/dt));%INIT
    NumberOfObs=0; % a running parametrs for the number of obstacles
    for ScanAngle = (AngleView-BeamWidth/2):dt:(AngleView+BeamWidth/2)
        SinScanAngle = sin(ScanAngle);
        CosScanAngle = cos(ScanAngle);
        for D = 0:DetectionDistance
            tempY = round(yCurrent + D.*SinScanAngle);
            tempX = round(xCurrent + D.*CosScanAngle);
            % Check the limits of terrain
            if tempY > LIMITS_1 || tempX > LIMITS_2 
                break
            elseif tempY <= 0 || tempX <=0
                break
            end %limits
            % find obstacles  and Prey in Terrain (the closest non-zeros at each alpha)
            CheckXY = Terrain(tempY,tempX);
            switch CheckXY
                case 0
                case 1 % obstacle
                    NumberOfObs=NumberOfObs+1;
                    ObsDistance(NumberOfObs) = D;
                    ObsAlpha(NumberOfObs) = ScanAngle-AngleView;
                    break % finds only the closest obstacle for each alpha
%                 case 2 % Prey -- NOT USED ___
%                     NumberOfPrey=NumberOfPrey+1;
%                     PreyDistance(NumberOfPrey) = D;
%                     PreyAlpha(NumberOfPrey) = ScanAngle-AngleView;
                otherwise
            end %switch CheckXY
        end %D
    end %ScanAngle
 
    % Outputs
    if NumberOfObs > 0
        ObsStruct.Distances = ObsDistance(1:NumberOfObs);
        ObsStruct.TargetAngle = ObsAlpha(1:NumberOfObs);
        IsObs = 1; % There is an obstacle
    else %if n=0  means it is free
        ObsStruct.Distances = [];
        ObsStruct.TargetAngle = [];
%             NearestDistance = [];
%             NearestAngle = [];
    end % if nOns
            ObsStruct.IsObs = IsObs;
end % BatFindEnvironment Function


