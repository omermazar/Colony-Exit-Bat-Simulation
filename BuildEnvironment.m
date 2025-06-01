function [ter,Xmax,Ymax,Zmax, Xmin, Ymin, Zmin, terrainPoints] = BuildEnvironment(TerrainParameters,xyResolution,EnvironmentType)

% BuildEnvironment MATLAB code for BatGUI1
%      BuildEnvironment, creates the terrain that the bat is flying in
%
%      Terrain = BuildEnvironment(TerrainParameters,SimulationParameters,EnvironmentType) returns 
%       the terrain the bat is flaying in;
%  Inputs:
%           TerrainParameters: a STRUCT with the following vaiables:
%                   Xmax, Ymax, Zmax - the size of the room [m]
%                   OBbsXWidth, OBbsYWidth - vectors with the sizes of
%                   rect. obsticles [m]
%                   OBbsDownLeftY, OBbsDownLeftX - Vectors with the position of the downleft corner
%           xyResolution: 
%                   the Resolution of the mesured room [m]  
%           EnvironmentType: a string with the following values
%                  'Free Space' - 10m room 
%                  'Empty Room'
%                   'Room 1' - a room with one obsticle
%                   'Room 2' - a room with two obsticles
% Outputs:
%           Terrain is a M x N matrix  
%                   M - Xmax/xyResolution , N - Ymax/xyResolution
%                   if Terrain(x,y) = 1 - there is an obsticle in (x,y)
%                   coordinates
%                   if Terrain(x,y) = 0- this coordinates is free
%                     Xmax, Ymax,Zmax- the size of the terrain

% Xmax = TerrainParameters.Xmax;
% Ymax = TerrainParameters.Ymax;
% Zmax = TerrainParameters.Zmax;
%modified by Johnathan 10.8.21
%%% Omer Nov2021

if isfield(TerrainParameters, 'TerrainMethod')
    TerrainMethod = TerrainParameters.TerrainMethod;
else
    TerrainMethod = 'byPoly';
end

Xmax = TerrainParameters.Xmax;
Ymax = TerrainParameters.Ymax;
Zmax = TerrainParameters.Zmax;
Xmin = 0; Ymin =0; Zmin =0;
        
        %% Walls
        shapeX = Xmax .* [0 0 1 1 0] / xyResolution ;
        shapeY = Ymax .* [0 1 1 0 0] / xyResolution  ;
        ter=[shapeX;shapeY];
        
switch EnvironmentType
    case 'Free Space'
        %in order to not confuse the program, the statment shouldnt be
        %empty, working as empty room for now
% %         shapeX = [0 0 0 0 0];
% %         shapeY = [0 0 0 0 0];
% %         ter=[ter ; shapeX/xyResolution ; shapeY/xyResolution ];
        shapeX = 100 .* [-1 -1 1 1 -1] / xyResolution ;
        shapeY = 100 .* [-1 1 1 -1 -1] / xyResolution  ;
        ter=[shapeX;shapeY];
    case 'Empty Room'
        %in order to not confuse the program, the statment shouldnt be
        %empty 
        shapeX = [0 0 0 0 0];
        shapeY = [0 0 0 0 0];
        ter=[ter ; shapeX/xyResolution ; shapeY/xyResolution ];
    case 'Room1' 
        %% Shapes
        n = 50; %Number of "leaves"
        SDleaf = 0.01; %SD of the shrub in metres 
        shapeX = [ 0 0 10 10 zeros(1,n) 0 ];
        shapeY = [ 0 10 10 0 zeros(1,n) 0 ];
        ter=[shapeX;shapeY];
        midX=5;
        midY=5;
        %Creating a bush
        shrub = randn(1,n)*SDleaf;
        shrubX = midX*ones(1,n) + shrub;
        shrub = randn(1,n)*SDleaf;
        shrubY = (midY-4)*ones(1,n) + linspace(3*SDleaf,8-3*SDleaf,n) + shrub;
        shapeX = [midX shrubX midX midX+0.5 midX+0.5 midX];
        shapeY = [midY-4 shrubY midY+4 midY+4 midY-4 midY-4];
        ter=[ter ; shapeX ; shapeY ]/xyResolution;
    case 'Room2'
        %% Shapes
        shapeX = [2 2 4 4 2];
        shapeY = [1 3 3 1 1];
        ter=[ter ; shapeX/xyResolution ; shapeY/xyResolution ];
        shapeX = [4 4 7 7 4];
        shapeY = [4 8 8 4 4];
        ter=[ter ; shapeX/xyResolution ; shapeY/xyResolution ];
    case 'Star'
        %in order to not confuse the program, the statment shouldnt be
        %empty 
        shapeX = [0 0 10 10 0 0 0 0 0 0 0 0 0];
        shapeY = [0 10 10 0 0 0 0 0 0 0 0 0 0];
        ter    = [shapeX/xyResolution;shapeY/xyResolution];
        midX   = 5;
        midY   = 5;
        shapeX = [midX midX+1 midX+2 midX+1.5 midX+2 midX+1 midX midX-1 midX-2 midX-1.5 midX-2 midX-1 midX];
        shapeY = [midY+2 midY+1 midY+1 midY midY-1 midY-1 midY-2 midY-1 midY-1 midY midY+1 midY+1 midY+2];
        ter    = [ter ; shapeX/xyResolution ; shapeY/xyResolution ];
   
    case 'semi-random obs'
         %% Shapes
        shapeX = [0 0 0 0 10 10 0 0 0 0 0 0 0 0 0 0 0];
        shapeY = [0 0 0 10 10 0 0 0 0 0 0 0 0 0 0 0 0];
        ter=[shapeX/xyResolution;shapeY/xyResolution];
        midX=5;
        midY=5;
        shapeX = [midX midX+0.5 midX+1.2 midX+1.7 midX+2.3 midX+1.2 midX+1.5 midX+1.1 midX midX-0.8 midX-0.9 midX-1.9 midX-1.2 midX-1.5 midX-2 midX-1.2 midX];
        shapeY = [midY+2.4 midY+2 midY+2.8 midY+1 midY+0.8 midY midY-1.1 midY-1.2 midY-2.1 midY-1.1 midY-1.3 midY-1 midY midY+1 midY+1.2 midY+1.3 midY+2.4];
        ter=[ter ; shapeX/xyResolution ; shapeY/xyResolution ];

    case 'Room-L'
         %% Shapes
         shapeLeftWallX = [Xmin Xmin Xmin+1 Xmin+1 Xmin];
         shapeLeftWallY = [Ymin Ymax-3 Ymax-3 Ymin Ymin];
         ter=[ter ; shapeLeftWallX/xyResolution ; shapeLeftWallY/xyResolution ];

         shapeRightUpX = [Xmax Xmax Xmin+3 Xmin+3 Xmax];
         shapeRightUpY = [Ymin+2 Ymax-3 Ymax-3 Ymin+2 Ymin+2];
         ter=[ter ; shapeRightUpX/xyResolution ; shapeRightUpY/xyResolution ];

         LeftDoorX = [Xmin+1 Xmin+1 Xmin+1.5 Xmin+1.5 Xmin+1];
         LeftDoorY = [Ymax-3 Ymax-3.5 Ymax-3.5  Ymax-3  Ymax-3];
         ter=[ter ; LeftDoorX/xyResolution ; LeftDoorY/xyResolution ];

         exitGap = 1;
         RightDoorX = [Xmin+2.5 Xmin+2.5 Xmin+3 Xmin+3 Xmin+2.5];
         RightDoorY = LeftDoorY;
         ter=[ter ; RightDoorX/xyResolution ; RightDoorY/xyResolution ];

    case 'Room-U'
        %% Shapes

        shapeMiddleWallX = [Xmax Xmax Xmax-8 Xmax-8 Xmax];
        shapeMiddleWallY = [Ymin+2 Ymax-2 Ymax-2 Ymin+2 Ymin+2];
        ter=[ter ; shapeMiddleWallX/xyResolution ; shapeMiddleWallY/xyResolution ];

        BottomDoorX = [Xmax-3 Xmax-3 Xmax-3.5 Xmax-3.5 Xmax-3];
        BottomDoorY = [Ymax-2 Ymax-1.5 Ymax-1.5  Ymax-2  Ymax-2];
        ter=[ter ; BottomDoorX/xyResolution ; BottomDoorY/xyResolution ];

        exitGap = 1;
        UptDoorX = BottomDoorX;
        UpDoorY = [Ymax Ymax-0.5 Ymax-0.5  Ymax  Ymax];
        ter=[ter ; UptDoorX/xyResolution ; UpDoorY/xyResolution ];
    
    case 'Room-L obs'
        %%
        %%% The Maze
        xMaze = [1.75 1.75 1.0  0.5  0.5 2.5 9.5 9.5 4 3.5 3.50 3.25 2.5  2.5  4    4   10  10 0 0   1.75] / xyResolution;
        yMaze = [9.3  8.75 8.75 8.0  2.5 0.5 0.5 3   3 3.5 8.00 8.75 8.75 9.3  9.3  3.5 3.5 0  0 9.3 9.3] / xyResolution;

        %%%% outer Walls
        shapeX = [-3 -3 12 12 -3*ones(1, numel(xMaze)-4)] / xyResolution ;
        shapeY = [-3 12 12 -3 -3*ones(1, numel(xMaze)-4)] / xyResolution ;
        ter=[shapeX;shapeY];
        ter = [ter; xMaze; yMaze];

        obsX = [1.5 1.5  2.75 2.75 1.5*ones(1, numel(xMaze)-4) ] / xyResolution;
        obsY = [5.5 6    6    5.5  5.5*ones(1, numel(xMaze)-4) ] / xyResolution  ;
        ter=[ter ; obsX ; obsY];
        
    case 'Simple Cave'
         %%
        %%% The Maze
        xMaze = [4.5  4.5   1.0  0.25  0.25 1    9.0  9.75 9.75 9    5.5  5.5 10  10 0  0  4.5 ] / xyResolution;
        yMaze = [10   9.75  9.75 9.0   1    0.25 0.25 1.0  9    9.75 9.75 10  10  0  0  10 10  ] / xyResolution;

        %%%% outer Walls
        shapeX = [-3 -3 12 12 -3*ones(1, numel(xMaze)-4)] / xyResolution ;
        shapeY = [-3 12 12 -3 -3*ones(1, numel(xMaze)-4)] / xyResolution ;
        ter=[shapeX;shapeY];
        ter = [ter; xMaze; yMaze];

        obsX = [4.5 4.5  5.5 5.5 4.5*ones(1, numel(xMaze)-4) ] / xyResolution;
        obsY = [5.5 6    6    5.5  5.5*ones(1, numel(xMaze)-4) ] / xyResolution  ;
        ter=[ter ; obsX ; obsY];
  

    otherwise
        %%
        error('unknown EnvironmentType');
end % switch EnvironmentType

%% Convert the Enironmemt to Points

%%% New Feb2023
if strcmp(TerrainMethod, 'byPoint')
    ds = TerrainParameters.TerrainGrid / xyResolution;
    terrainPoints =  polyToPoints(ter, ds);
    %%%% New June 2024
    % remove points from the final Exit
    switch EnvironmentType
        case 'Room-L obs'
            ixUpper = terrainPoints(1,:) >= 0 & terrainPoints(1,:) <= 400 & terrainPoints(2,:) >= 900;
            terrainPoints(:,ixUpper) = nan;
    end

else
   terrainPoints  = [];
end

%%% Old Definitions

%   case 'Room-L'
%         T_W      = 4;
%         shapeX   = [Xmin+1 Xmin  Xmin  Xmax  Xmax T_W  T_W   T_W-1  T_W-1   T_W+0.25   T_W+0.25   Xmax+0.25  Xmax+0.25 Xmin-0.25 Xmin-0.25 Xmin+1 Xmin+1];
%         shapeY   = [Ymax   Ymax  Ymin  Ymin  T_W  T_W  Ymax  Ymax   Ymax+0.25  Ymax+0.25  T_W+0.25    T_W+0.25  Ymin-0.25 Ymin-0.25 Ymax+0.25 Ymax+0.25 Ymax];
%         %         xf = [ Xmin+1   Xmin+1    T_W-1   T_W-1  Xmin+1];
%         %         yf = [ Ymax     Ymax+0.25    Ymax+0.25  Ymax   Ymax ];
%         zeroFill = zeros(2, numel(shapeX) - size(ter,2) );
%         ter      = [[ter, zeroFill] ; shapeX/xyResolution ; shapeY/xyResolution ];
% 
%     case 'Room-U'
%         T_W      = 4;
%         shapeX   = [Xmax Xmax   Xmin-1 Xmin-1 Xmax+1 Xmax+1 T_W+1 T_W+1      Xmax       Xmax     T_W       T_W Xmax Xmax Xmin Xmin Xmax];
%         shapeY   = [Ymax Ymax+1 Ymax+1 Ymin-1 Ymin-1 T_W+1  T_W+1 Ymax-T_W-1 Ymax-T_W-1 Ymax-T_W Ymax-T_W  T_W T_W  Ymin Ymin Ymax Ymax];
%         zeroFill = zeros(2, numel(shapeX) - size(ter,2) );
%         ter      = [[ter, zeroFill] ; shapeX/xyResolution ; shapeY/xyResolution ];