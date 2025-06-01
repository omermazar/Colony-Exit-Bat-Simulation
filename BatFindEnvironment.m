
function[ ObsStruct, PreyStruct, PreyDist, PreyAngle, Bat2PreyRelativeAngle] = ...
        BatFindEnvironment(xCurrent,yCurrent,AngleView,BeamWidth, ...
        DetectionDistance, Terrain, PreysPositions)
    
% This funtion will find if there are an obsticles or preys in the detection range...
% and return the range and angle from the current point
% Outputs:
    % IsObs = 1 if there was an obstacle Found
    % ObsStruct - structure of the obstacle found:
    %       ObsStruct.Distance is the a vector of Distances between the bat and the
    %       obstacles found
    %       ObsStruct.Alpha are the angles between the bat flight direction and the obstacle
    %       ObsStruct.IsObs -1 if there was an obstacle Found 
    % NumOfPreysDetected = 1the number of the prys that eat detected
    % DetectedPreys - which Preys were detected
    % PreyStruct - structure of the prey found:
    %       PreyStruct.Prey(n) - the n-th prey is struct with the fields: 
    %            struct('IsDetected',, 'Distance',[], 'Alpha',[], 'RelativeDirection',[])
    %       PreyStruct.Alpha is the angle between the bat flight direction and the obstacle
    %       PreyStruct.IsPrey - 1 if there was a prey
    % PreyDist - vector of distances from the bat to all preys
    % PreyAngle - vector of angles from the bat to all preys
    
% Inputs: 
    % xCurrent,yCurrent, - position and direction of the bat
    % AngleView - The Direction of the Bat Head (for now equals to his direction of
    % movement)
    % BeamWidth, DetectionDistance - The detection Sonar Beam-Width and
    % distance in the terrain coordinations system (no-units)
    % Terrain = the environment to check
    % PreysPositions = [xPrey,yPrey,PreyTeta, PreyVelictu] X Number Of Preys 
    
    %%% INIT %%%
    IsObs=0; IsPreyDetected=0; % Init
    MatSize = size(PreysPositions);
    TotalNumberOfPreys = MatSize(1);
    NumOfPreysDetected = 0;
    DetectedPreys = zeros(1,TotalNumberOfPreys); 
    dt=0.1; % the step size for the angle in radians
%     PreyStruct.IsPrey = 0;
    PreyStruct.Prey(TotalNumberOfPreys) = ...
        struct('IsDetected',0, 'Distance',[], 'Alpha',[], 'RelativeDirection',[], 'Bat2PreyRelativeAngle',[]);
%     PreyStruct.Distance = [];
%     PreyStruct.Alpha = [];
%     PreyStruct.RelativeDirection = [];
    
    %%%% Find the Prey %%%
   
    if TotalNumberOfPreys > 0
        PreyDist = zeros(1,TotalNumberOfPreys);
        PreyAngle =  zeros(1,TotalNumberOfPreys);
        Bat2PreyRelativeAngle = zeros(1,TotalNumberOfPreys);
        for PreyNum = 1:TotalNumberOfPreys
            
            xPrey = PreysPositions(PreyNum,1);
            yPrey = PreysPositions(PreyNum,2);
            PreyTeta= PreysPositions(PreyNum,3);
            [PreyDist(PreyNum), PreyAngle(PreyNum), Bat2PreyRelativeAngle(PreyNum)] = ...
                FindDistAngle(xCurrent, yCurrent, AngleView, xPrey, yPrey, PreyTeta);

            IsPreyInRange = (PreyDist(PreyNum) <= DetectionDistance);
            if IsPreyInRange
                if ( PreyAngle(PreyNum) <  BeamWidth/2 &&  PreyAngle(PreyNum) > -BeamWidth/2 )
                    IsPreyDetected =1; % FLAG if the Bat Detects the Prey
                    %%% Update the outopt struct %%%
                    PreyStruct.Prey(PreyNum).IsDetected = 1;
                    PreyStruct.Prey(PreyNum).Distance = PreyDist(PreyNum);
                    PreyStruct.Prey(PreyNum).Alpha = PreyAngle(PreyNum);
                    PreyStruct.Prey(PreyNum).Bat2PreyRelativeAngle = Bat2PreyRelativeAngle(PreyNum);
                    %                         PreyStruct.Prey(PreyNum).RelativeDirection = PreyFlightDirection - AngleView;
                    NumOfPreysDetected = NumOfPreysDetected+1;
                    DetectedPreys(NumOfPreysDetected) = PreyNum;
                end
            end  %IsPreyDist
    % %         end %  if PreyDist <  
    %         Terrain(round(yPrey) , round(xPrey)) = 2; % if there is a prey the terrain is 2 at this pos
        end %for PreyNum = 1:TotalNumberOfPreys
        PreyStruct.DetectedPreys= nonzeros(DetectedPreys)';
        PreyStruct.NumOfPreyDetected = NumOfPreysDetected;
    end % if TotalNumberOfPreys > 0
    
    
    %%%% find Obstacles %%%
    
    LIMITS = size(Terrain);
    ObsDistance=zeros(1,ceil(DetectionDistance/dt)); % INIT
    ObsAlpha=zeros(1,ceil(BeamWidth/dt));%INIT
    PreyDistance=zeros(1,ceil(DetectionDistance/dt)); % INIT
    PreyAlpha=zeros(1,ceil(BeamWidth/dt));%INIT
    NumberOfObs=0; % a running parametrs for the number of obstacles
    NumberOfPrey=0;  % a running parametrs for the number of Preys
    for ScanAngle = (AngleView-BeamWidth/2):dt:(AngleView+BeamWidth/2) 
        for D = 0:DetectionDistance
            tempY = round(yCurrent+D.*sin(ScanAngle));
            tempX = round(xCurrent+D.*cos(ScanAngle));
            % Check the limits of terrain
            if (tempY > LIMITS(1)) || (tempX > LIMITS(2)) 
                break
            elseif (tempY <= 0) || (tempX <=0)
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
                case 2 % Prey -- NOT USED ___
                    NumberOfPrey=NumberOfPrey+1;
                    PreyDistance(NumberOfPrey) = D;
                    PreyAlpha(NumberOfPrey) = ScanAngle-AngleView;
                otherwise
            end %switch CheckXY
        end %D
    end %ScanAngle
 
    % Outputs
    if NumberOfObs > 0
        ObsStruct.Distance = ObsDistance(1:NumberOfObs);
        ObsStruct.Alpha = ObsAlpha(1:NumberOfObs);
        IsObs = 1; % There is an obstacle
%             ObsDistaneFromLinear = ObsDistance.*abs(tan(ObsAlpha)); % for calculates if we neend to mauver 
%             [NearestDistance,nn] = min(min(ObsDistaneFromLinear, ObsDistance));
%             NearestAngle = ObsAlpha(nn);
    else %if n=0  means it is free
        ObsStruct.Distance = [];
        ObsStruct.Alpha = [];
%             NearestDistance = [];
%             NearestAngle = [];
    end % if nOns

            ObsStruct.IsObs = IsObs;
%             PreyStruct.IsPrey = IsPreyDetected;
end % BatFindEnvironment Function

% % 
% % function [R] = findHipDist(xCurrent,yCurrent,PreyDist,PreyAngle,xPrey,yPrey)
% %     XPreyHip =  xCurrent + PreyDist*cos(PreyAngle);
% %     YPreyHip =  yCurrent + PreyDist*sin(PreyAngle);
% %     R = sqrt((XPreyHip-xPrey)^2 + (YPreyHip-yPrey)^2 );
% % end
% % 
% %             
