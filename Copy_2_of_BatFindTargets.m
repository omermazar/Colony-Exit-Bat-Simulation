
function [TargetStruct] = ...
        BatFindTargets(xCurrent,yCurrent,AngleView,BeamWidth, ...
        DetectionDistance, Terrain, TargetPositions, TargetType)
    
% This funtion will find if there are an obstacles or preys in the detection range...
% and return the range and angle from the current point
% Outputs:
    % DetectedPreys - which Preys were detected
    % TargetStruct - structure of the targets found:
    %       TargetsID - a vector of the Preys ID (numbers)in the area
    %       NumOfTagets: TotalNumber of preys in the area
    %       Distances - a vector of the disace fro the bat to each prey 
    %       TargetAngle: vector of the angles
    %       Bat2PreyRelativeAngle - vector of the relative direction
    %               between the bat movemnet and the target movement
    %       IsInBeam - vector of zeros (no) or ones(yes)

    
% Inputs: 
    % xCurrent,yCurrent, - position and direction of the bat
    % AngleView - The Direction of the Bat Head (for now equals to his direction of
    % movement)
    % BeamWidth, DetectionDistance - The detection Sonar Beam-Width and
    % distance in the terrain coordinations system (no-units)
    % Terrain = the environment to check
    % PreysPositions = [xPrey,yPrey,PreyTeta, PreyVelictu] X Number Of Preys 
    % TargetType  = 'Prey' or 'OtherBats'
    
    %%% INIT %%%
   
    MatSize = size(TargetPositions);
    NumOfTargets = MatSize(1);
    
% % %     switch TargetType
% % %         case 'Prey'
% % %             NumOfTargets = MatSize(1);
% % %         case 'OtherBats'
% % %             NumOfTargets = MatSize(1) - 1;
% % %     end % switch TargetType

    TargetStruct = ...
        struct('NumOfTargets',NumOfTargets, 'TargetsID',[], 'Distances',[], 'TargetAngle',[], 'Bat2TargetRelativeAngle',[], 'IsInBeam',[]);
    
    %%%% Find the Prey %%%
   
    if NumOfTargets > 0
        TargetStruct.Distances = zeros(1,NumOfTargets);
        TargetStruct.TargetAngle =  zeros(1,NumOfTargets);
        TargetStruct.Bat2TargetRelativeAngle = zeros(1,NumOfTargets);
        TargetStruct.IsInBeam = zeros(1,NumOfTargets);
        %         TargetStruct.TargetsID = 1:NumOfTargets;
        for TargetNum = 1:NumOfTargets
            %             TargetNum = TargetPositions(TargetNum,1);
            TargetStruct.TargetsID(1,TargetNum) = TargetPositions(TargetNum,1);
            xPrey = TargetPositions(TargetNum,2);
            yPrey = TargetPositions(TargetNum,3);
            PreyTeta= TargetPositions(TargetNum,4);
            [TargetStruct.Distances(TargetNum), TargetStruct.TargetAngle(TargetNum),...
                TargetStruct.Bat2TargetRelativeAngle(TargetNum)] = ...
                FindDistAngle(xCurrent, yCurrent, AngleView, xPrey, yPrey, PreyTeta);
            %%% Check if in Beam
            if (TargetStruct.Distances <= DetectionDistance)
                if ( TargetStruct.TargetAngle(TargetNum) <  BeamWidth/2 )...
                        && ( TargetStruct.TargetAngle(TargetNum) > -BeamWidth/2 )
                    TargetStruct.IsInBeam(TargetNum) = 1;
                end % if ( TargetAngle(PreyNum) <  BeamWidth/2 &&  TargetAngle(PreyNum) > -BeamWidth/2 )
            end %if (TargetDist(PreyNum) <= DetectionDistance)
            
            
        end % for TargetNum = 1:NumOfTargets
            
        
    end %  if NumOfTargets > 0

