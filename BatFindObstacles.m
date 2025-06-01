function [ObsStruct] = BatFindObstacles(xBat,yBat,tetaBat, BatBeamWidth, DetectionRange, terrain, polyTerrain, AllParams  )

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
    % polyTerrain: 
    %               for k = 1:2:size(Terrain,1)-1
                        %     polyTerrain{k} = polyshape( Terrain(k,:),Terrain(k+1,:));
                        % end

%% Input
if isfield(AllParams.TerrainParams, 'TerrainMethod')
    TerrainMethod = AllParams.TerrainParams.TerrainMethod;
else
    TerrainMethod = 'byPoly';
end

%%% INIT %%%
%%% Vladimir
% Omer - Nov2021
% Major Bug FIX 14Jun22

ObsStruct.Distances            = [];
ObsStruct.TargetAngle          = [];
ObsStruct.xOBSTC               = [];
ObsStruct.yOBSTC               = [];
ObsStruct.IsObs                = 0;
ObsStruct.DetectetedTargetsVec = [];


%% Object detect loop
% if free space- no needto calculate
switch TerrainMethod
    %% Find Obasltcles by intersections between the bema-lines and the polyons of the environment 
    case 'byPoly'

    %%% Resolution - by how much to devide the BatBeamWidth
    resolution = 2*pi/180; % 1.5 [deg]
    sz = size(terrain,1);
    nPoints =  DetectionRange /100 / 0.01 ; % resolutin of 1cm in detection

    if ~strcmp(AllParams.TerrainParams.Type, 'Free Space')
        for t=tetaBat-BatBeamWidth/2:+resolution:tetaBat+BatBeamWidth/2
            % Search Line
            srchlineX = [xBat, xBat + DetectionRange*cos(t)]; % linspace(xBat, xBat+DetectionRange*cos(t), nPoints);
            srchlineY = [yBat, yBat + DetectionRange*sin(t)]; % linspace(yBat, yBat+DetectionRange*sin(t), nPoints);
            %lntersecId=0

            % Intersection border
            if ~strcmp(AllParams.TerrainParams.Type, 'Room-L obs')
                %%% in 'Room-L obs' the bat is inside the cave, no need to
                %%% check room-borders
                %             [in, on] = inpolygon(srchlineX, srchlineY, terrain(1,:), terrain(2,:));
                %             poly2 = polyshape(terrain(1,:), terrain(2,:));
                [in, out] = intersect(polyTerrain{1}, [srchlineX; srchlineY]');

                %             lntersecId = find(in==0,1,'first');
                if any(out) % lntersecId>0
                    ObsStruct.xOBSTC = [ObsStruct.xOBSTC out(1,1)]; %[ObsStruct.xOBSTC srchlineX(lntersecId)];
                    ObsStruct.yOBSTC = [ObsStruct.yOBSTC out(1,2)];% [ObsStruct.yOBSTC srchlineY(lntersecId)];
                    try
                        %                     ObsStruct.Distances = [ObsStruct.Distances pdist([xBat,yBat;srchlineX(lntersecId),srchlineY(lntersecId)],'euclidean')];
                        ObsStruct.Distances = [ObsStruct.Distances pdist([xBat,yBat; out(1,1),out(1,2)],'euclidean')];
                    catch
                        warning('oooopppsssss: pdist on BatFindObstacles')
                    end % try
                    ObsStruct.TargetAngle = [ObsStruct.TargetAngle, t-tetaBat];
                    ObsStruct.IsObs=1;
                end
            end     %   if ~strcmp(AllParams.TerrainParams.Type, 'Room-L obs')


            % Intersection objects
            for i=3:+2:sz-1 %for every object loop
                %             [in, on] = inpolygon(srchlineX, srchlineY, terrain(i,:), terrain(i+1,:));
                %             poly2 = polyshape(terrain(i,:), terrain(i+1,:));
                [in, out] = intersect(polyTerrain{i}, [srchlineX; srchlineY]');
                addObjFlg = false;
                %             lntersecId = find(in,1,'first');
                %             if lntersecId>0

                if any(in) % lntersecId>0
                    currTeta = t-tetaBat;
                    currDist = pdist([xBat,yBat; in(1,1),in(1,2)],'euclidean');
                    [~, ixTeta]   = ismembertol(currTeta, ObsStruct.TargetAngle, resolution*0.25, 'DataScale', 1);
                    %%% add to detection only if the current Distance is smaller
                    %%% than the Distance in the detection vector
                    if ixTeta == 0
                        % New teta - Add
                        addObjFlg = true;
                        ObsStruct.Distances   = [ObsStruct.Distances, currDist];
                        ObsStruct.TargetAngle = [ObsStruct.TargetAngle, currTeta];
                        ObsStruct.xOBSTC      = [ObsStruct.xOBSTC in(1,1)]; %[ObsStruct.xOBSTC srchlineX(lntersecId)];
                        ObsStruct.yOBSTC      = [ObsStruct.yOBSTC in(1,2)];% [ObsStruct.yOBSTC srchlineY(lntersecId)];
                        ObsStruct.IsObs       = 1;

                    elseif currDist < min(ObsStruct.Distances(ixTeta))
                        % Teta exists and Distance is closer- Replace
                        addObjFlg = true;
                        ObsStruct.Distances(ixTeta)   = currDist;
                        ObsStruct.TargetAngle(ixTeta) = currTeta;
                        ObsStruct.xOBSTC(ixTeta)      = in(1,1); %[ObsStruct.xOBSTC srchlineX(lntersecId)];
                        ObsStruct.yOBSTC(ixTeta)      = in(1,2);% [ObsStruct.yOBSTC srchlineY(lntersecId)];
                        ObsStruct.IsObs               = 1;

                    else %
                        % Teta exists and Distance is not closer- Do Nothing
                        addObjFlg = false;
                    end % if ~ismembertol()


                end % if any(in)
            end % for i = 3: ..

        end % for t = teta: ...
    end % if ~strcmp(AllParams.TerrainParams.Type, 'Free Space')
    
    %% Find Obasltcles by the points of the environment 
    case 'byPoint'
        %% Detectable points inside the beam
        terrainPoints = AllParams.TerrainParams.terrainPoints;
        in = sqrt((terrainPoints(1,:) - xBat).^2 + (terrainPoints(2,:)-yBat).^2) <= DetectionRange;
        % %% points in range
        pointsInBeam.xin = terrainPoints(1,in); 
        pointsInBeam.yin = terrainPoints(2,in);
        xin = pointsInBeam.xin; yin = pointsInBeam.yin;
        % Distances and Angles
        pointsInBeam.distances = sqrt((xBat-xin).^2 + (yBat-yin).^2);
        pointsInBeam.tetas     = wrapToPi( atan2((yin-yBat), (xin-xBat)) - tetaBat );

        %%%% Remove points that are hidden by another wall
        minDist = AllParams.TerrainParams.TerrainGrid / AllParams.SimParams.xyResolution * 1.1;
   
        if ~isempty(pointsInBeam)
            ixPoints = filterFarSide(pointsInBeam, minDist,0);
            
            %%% find the pionts that are within the actual beam
            inBeam = pointsInBeam.tetas(ixPoints) <= (BatBeamWidth/2) & ...
                pointsInBeam.tetas(ixPoints) >= -(BatBeamWidth/2);
            ixPoints = ixPoints(inBeam);    
        else
            ixPoints = [];
        end % pointsInBeam
 
        %%% Output
        ObsStruct.TargetAngle = pointsInBeam.tetas(ixPoints);
        ObsStruct.Distances   = pointsInBeam.distances(ixPoints);
        ObsStruct.xOBSTC      = pointsInBeam.xin(ixPoints);
        ObsStruct.yOBSTC      = pointsInBeam.yin(ixPoints);
        ObsStruct.IsObs       = any(ixPoints);
end% swtich terrainMethod

end % function

