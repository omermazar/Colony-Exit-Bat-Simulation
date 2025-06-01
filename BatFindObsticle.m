
function[IsFree,Distance,Alpha,NearestDistance,NearestAngle] = ...
        BatFindEnvironment(xCurrent,yCurrent,AngleView,BeamWidth, ...
        DetectionDistance, Terrain)
    
% This funtion will find if ther is an obsticle in the detection range...
% and return the range and angle from the current point
% Outputs:
    % Isfree = 1 if there was no obstacle Found
    % Distance is the Distance between the bat and the obstacle
    % Alpha is the angle between the bat flight direction and the obstacle
    % NearestDitance is the Ditnce to the nearst obstacle from the bat or frrom
    % the bat direction flight
    % NearestAngle is the angle between the bat and the nearest obstacle  
% Inputs: 
    % xCurrent,yCurrent,AngleView - position and direction of the bat
    % BeamWidth, DetectionDistance - THe detection Sonar Beam-Width and
    % distance
    % global Terrain;
    
    dt=0.1; % the step size for the angle in radians
    Distance=zeros(1,ceil(BeamWidth/dt)); % INIT
    Alpha=zeros(1,ceil(BeamWidth/dt));%INIT
    n=0;  
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
            % find obstacles in Terrain (the closest non-zeros at each alpha)
           if Terrain(tempY,tempX) ~= 0
              n=n+1;
              Distance(n) = D;
              Alpha(n) = ScanAngle-AngleView;
              break % finds only the closest obstacle for each alpha
           else
                
           end %if Terrain()
        end %D
    end %ScanAngle
    % Outputs
    if n > 0
        Distance = Distance(1:n);
        Alpha = Alpha(1:n);
        IsFree = 0; % There is an obstacle
        ObsDistaneFromLinear = Distance.*abs(tan(Alpha)); % for calculates if we neend to mauver 
        [NearestDistance,nn] = min(min(ObsDistaneFromLinear, Distance));
        NearestAngle = Alpha(nn);
    else %if n=0  means it is free
        Distance = [];
        Alpha = [];
        IsFree = 1; % There is NO obstacle
        NearestDistance = [];
        NearestAngle = [];
    end % if n
  
end % FindObstacle Function

