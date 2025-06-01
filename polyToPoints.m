
function [terrainPoints] =  polyToPoints(Terrain, ds) 
% [Terrain, TerrainParams.Xmax, TerrainParams.Ymax, TerrainParams.Zmax,...
%     TerrainParams.Xmin, TerrainParams.Ymin, TerrainParams.Zmin] = ...
%     BuildEnvironment(AllParams.TerrainParams, AllParams.SimParams.xyResolution, AllParams.TerrainParams.Type)

% ds = 10;
xEdge =[]; yEdge = [];

% each doublle row in trrain is a polygon
for kRow = 1: 2: (size(Terrain,1)-1)
    % find each line angles 
    for nPoint = 1:(size(Terrain,2)-1)
        xStart = Terrain(kRow,nPoint);
        xEnd   = Terrain(kRow,nPoint+1);
        yStart = Terrain(kRow+1,nPoint);
        yEnd   = Terrain(kRow+1,nPoint+1);

        lineAngle = atan2(yEnd - yStart, xEnd - xStart);
        nPoints   = floor(sqrt((xEnd-xStart).^2 + (yEnd- yStart).^2) / ds);
        
        xL = xStart + cumsum(ds*cos(lineAngle)*ones(1,nPoints));
        yL = yStart + cumsum(ds*sin(lineAngle)*ones(1,nPoints));
        
        xEdge = [xEdge, xL, xEnd];
        yEdge = [yEdge, yL, yEnd];
    end % for nPoint
end % for kPoly

terrainPoints = [xEdge; yEdge];

% remove double points (corners)
ix = diff(terrainPoints(1,:)).^2 + diff( terrainPoints(2,:).^2) == 0 ;
terrainPoints(:,ix) = [];

