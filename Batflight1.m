function [Bat_Position_Vec] = Batflight()
% THIS Function simulates and plots a bat flght in a room 
% Types of movement:
% P2P (Point 2 Point) - Constant Velocity, Constant angle (teta)
% Circles - Constant velocity and angular velocity (dteta)
% Eights - 
% RandomFlight (dvelocity, dteta are gausian variables)

% ObsticlesAVoidance - when the bat approches obsticle it manuvers to avoid
% the hit
% first we simulate 2d flight

%%%%%%%%%%%%%
% Simulation Paramaters
xyResolution = 0.01; % 0.01=1cm
SimulationTime = 20; % Sec
SampleTime = 0.001; % 0.001= 1miliSec
NumOfSaples = SimulationTime/SampleTime; % The total number of samples
TimeVec= [0:SampleTime:SimulationTime];
AvoidObsticleFlag= true; % to manuver away from absticles
AnimationFlag = true; % Animate Bat Movement in time
FlightTypeFlag = 'circle'; % 'random' or 'linear' or 'circle'
%%%%%%%%%%%%%

%%%%%%%%%%%%%%
% The TERRAIN Parameters- closed room 10X10X10m with Rectanguler Obsticle
% if Terrain(x,y)= 0 - there is no obsticle at x,y
% if Terrain(x,y)= 1 - there is an obsticle at x,y
Xmax = 15; Ymax = 15; Zmax = 10; % meters
[Xcor, Ycor] = meshgrid(0:xyResolution:Xmax,0:xyResolution:Xmax);
Terrain = zeros( ceil(Xmax/xyResolution)+1 , ceil(Ymax/xyResolution)+1 );
Terrain(1,:)=1;  Terrain(end,:)=1;   Terrain(:,1)=1;   Terrain(:,end)=1; % Room Walls
LIMITS = size(Terrain);

% Rectangle Obsticle
XWidth = 0.5; YWidth = 1; % meters
DownLeftY = 2.5;DownLeftX = 3;% position of left buttom corner (meters)
xRect= [ ceil(DownLeftX/xyResolution) : (ceil(DownLeftX/xyResolution)...
+ ceil(XWidth/xyResolution)) ];
yRect= [ ceil(DownLeftY/xyResolution) : (ceil(DownLeftY/xyResolution)...
+ ceil(YWidth/xyResolution)) ];
Terrain(yRect,xRect) = 1;


% PLOT The Terrain 
close all
figure (1)
contour(Xcor,Ycor,Terrain)
xlabel('x [m]'); ylabel('y [m]');
title('The Room')
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
% BAT FLIGHT - Basic Parameters
MaxVelocity = 12 ; % m/s
NominalVelocity = MaxVelocity*0.25; % m/s 
MaxAccel = 20; % [m/sec^2] ~2g 
NominalAccel = MaxAccel*0.5; % m/s^2
MaxOmega = MaxAccel/NominalVelocity; % Maximum Radian Velocity [rad/sec]
% Radius of circular = Velocity/Omega
% R = NominalVelocity/MaxOmega % ~0.1m the Radius of the bat Manuver
Omega = NominalAccel/NominalVelocity;% [rad/sec] 
MinDistanceAllowedM = 0.3; % [m] The minimum range the bat can be near Obsticle

% BAT SONAR BASIC
BatBeamWidth = pi/6; % radians
BatDetectionDistance = 2; % meters
BatPRF = 25; % [1/sec] Trasmits per seconds
%%%%%%%%%%%%%%%%%%%%%%

% STARTING Point
StartX0 = 0.5; StartY0= 0.5; % meters
Teta0 = pi/3; % Direction in radian 
Accel0 = 0; % No accelaration
Dteta0 = 0; % Direct Line

%%%%%%%%%%
% Position, Velocity and Time in Simulation Index (no units)
X0 = StartX0/xyResolution;
Y0 = StartY0/xyResolution;
V0 = NominalVelocity*SampleTime/xyResolution;
DTETA = Omega*SampleTime;
MaxdTeta = MaxOmega*SampleTime;
DetectionDistance = BatDetectionDistance / xyResolution;
RegTimeToSonar = 1/BatPRF/SampleTime; % regular- in free terrain
ManuverTimeToSonar = 1/BatPRF/SampleTime/10; % PRF when Manuvering
MinDistanceAllowed = MinDistanceAllowedM/xyResolution;
% Initial V0,Xbati,Ybati
% (xBati, yBati) - The position of the Bat in the simulation
Velocity = V0;
xBati= zeros(1,NumOfSaples); 
yBati= zeros(1,NumOfSaples); % Initial
xBati(1) = X0; 
yBati(1) = Y0;
%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%
%%%%%%%%%%
% THE MAIN FUNCTION

Teta = Teta0*ones(1,NumOfSaples); % Initial Linear movement 
DTetaManuever = 0; % Initial
if (AvoidObsticleFlag)
    tTimeToSonar = 0; % Time to next search, at the start- SEARCH
    % Help Vectors for debuggind
        IsFreeVec = ones(1,NumOfSaples);
        xFind = [];
        yFind = [];     
       ManueverFlag = zeros(1,NumOfSaples);
%         ObsRangeVec = zeros(1,NumOfSaples);
%         ObsAngleVec = zeros(1,NumOfSaples);
    for nTime = [1:NumOfSaples]     %nTime- the time vector insample
        if tTimeToSonar == 0
            [IsFREE,ObsRange,ObsAngle] = FindObsticle(xBati(nTime),yBati(nTime),Teta(nTime)); % % DetectionDistance,BatBeamWidth);
            % IsFREE= 1 if there is no obsticle
            % (ObsRAnge,ObsAngle) - vector of the positions of the obstcles
            % .. range and angle from the bat
            if (~IsFREE)
                 % The Position of the Obsicles found
                IsFreeVec(nTime) = IsFREE;
                xFind = [xFind, xBati(nTime)+ObsRange.*cos(Teta(nTime)+ObsAngle)];
                yFind = [yFind, yBati(nTime)+ObsRange.*sin(Teta(nTime)+ObsAngle)];
                 % find the nearest obticle to avoid
                % nearest obs
                ObsDistaneFromLinear = ObsRange.*abs(tan(ObsAngle)); % for calculates if we neend to mauver 
                [ObsToAvoidDistance,nn] = min(min(ObsDistaneFromLinear, ObsRange));
                ObsToAvoidAngle = ObsAngle(nn);
                if ObsToAvoidDistance <= MinDistanceAllowed
                    ManueverFlag(nTime:nTime+ManuverTimeToSonar-1) = 1;  % if manuver is needed- keep manuver until the next pulse
                    DTetaManuever = -sign(ObsToAvoidAngle)*MaxdTeta;
                    tTimeToSonar = ManuverTimeToSonar; % while manuver allways check Distace
                else % if ObsToAvoidDistance
                    ManueverFlag(nTime) = 0;
                    DTetaManuever = 0;
                    tTimeToSonar = RegTimeToSonar; 
                end % if ObsToAvoidDistance
            end % if IsFree
        else % no search- continue to flight
%             ManueverFlag(nTime) = 0;
%             DTetaManuever =0;
            tTimeToSonar = tTimeToSonar -1;
        end % if tTimeToSonar
        
        % Movement of the BAT
            [xBati(nTime+1),yBati(nTime+1),Teta(nTime+1)] =  ...
             BatMovemet(xBati(nTime),yBati(nTime), Teta(nTime), FlightTypeFlag,...
             ManueverFlag(nTime), DTetaManuever, DTETA, V0); 

    end  % for nTime - Time     
elseif (~AvoidObsticleFlag)% AvoidObsticle  - open Space movement
    Teta = Teta0; 
    ManueverFlag = 0;
    DTetaManuever = 0;
    for nTime = [1:NumOfSaples]     %nTime- the time
            [xBati(nTime+1),yBati(nTime+1),Teta(nTime+1)] =  ...
             BatMovemet(xBati(nTime),yBati(nTime), Teta(nTime),FlightTypeFlag, ManueverFlag,...
             DTetaManuever,DTETA, V0); 
        end % for nTime   
end % if Avoid Flag


%%%%%%%%%
%%%%%%%%%
function[IsFree,Distance,Alpha] = FindObsticle(xCurrent,yCurrent,AngleView)
% This funtion will find if ther is an pbsticle in the detection range...
% and return the range and angle from the current point
% Isfree = 1 if there was no Obsticle Found
% Distance is the Distance between the bat and the obsticle
% Alph is the angle betwin the bat flight direction and the obsticle
    Distance=[];
    Alpha=[];
    for ScanAngle = (AngleView-BatBeamWidth/2):0.1:(AngleView+BatBeamWidth/2)
        for D = [0:DetectionDistance]
            tempY = round(yCurrent+D.*sin(ScanAngle));
            tempX = round(xCurrent+D.*cos(ScanAngle));
            % Check the limits of terrain
            if (tempY > LIMITS(1)) || (tempX > LIMITS(2)) 
                break
            elseif (tempY <= 0) || (tempX <=0)
                break
            end %limits
            % find Obsticles in Terrain (the closest non-zeros at each alpha)
           if Terrain(tempY,tempX) ~= 0
              Distance = [Distance, D];
              Alpha = [Alpha, ScanAngle-AngleView];
              break % finds only the closest obsticle for each alpha
           else
                
           end %if Terrain()
        end %D
    end %ScanAngle
    IsFree = isempty(Distance);
  
end % FindObsticle Function

%%%%%%%%%
%%%%%%%%%
function[xBNext,yBNext,TetaBNext] =  BatMovemet(xB,yB, TetaB, FlightType, ManFlag, DTetaMan, DTetaBat, Vel); 
% This funtion will calculate tne new position(x,y,teta) of the  bat
% output arguments: 
% xBNext -the x coordinate ,  yBNext -the y coordinate , 
% TetaBNext- the angle of the bat
% Inptus: 
% xB,yB, TetaB- current position
% Vel = teh velocity of the bat
% FlightType = 'linear', 'lircle', random'
% ManFlag, DTetaMan - flag weather to nmanuver, and teh differncial angle
% DTetaBat =- the diifferntial angle of the bat ' from its accelartion
 KdTeta = 5;  % 1/SampleTime; 
    if ManFlag == 0 % no need to manuver
        switch FlightType
            case 'linear'
                dTeta = 0; 
            case 'circle'
                dTeta =  DTetaBat*KdTeta;
            case 'random'
               dTeta = DTetaBat*randn*KdTeta*2; 
            otherwise
                disp('other value')
        end % switch FlightType
    elseif ManFlag == 1
        dTeta = DTetaMan*10;
    end; %ManFlag
    TetaBNext = TetaB+ dTeta;;
    xBNext = xB + Vel.*cos(TetaBNext) ;
    yBNext = yB + Vel.*sin(TetaBNext) ;
end% BatMovement Function


% Final DATA
xBatPos = xBati*xyResolution;
yBatPos = yBati*xyResolution;
if (AvoidObsticleFlag)
    xObsFinds = xFind*xyResolution;
    yObsFinds = yFind*xyResolution;
end % AvoidObs

% PLOTS
hold on
title('BAT Movement #1')
plot(xBatPos,yBatPos,'ok')
if (AvoidObsticleFlag)
    plot(xBatPos(find(~IsFreeVec)),yBatPos(find(~IsFreeVec)),'g*') % the postion when finfdingobsticle
    plot(xObsFinds,yObsFinds,'+r')
end % AvoidObs


% animation
 f1 = figure(2);
if AnimationFlag
    subplot(2,1,1)
    title('Bat Movement Animation')
    h = animatedline('Color','k'); % subplot(2,1,1)
    hprevoius = animatedline('Color','c'); % Show prev location
    axis([min(xBatPos),max(xBatPos),min(yBatPos),max(yBatPos)]);
    subplot(2,1,2)
    title('Bat Current Direction Animation')
    axis([0,SimulationTime,min(Teta)/pi,max(Teta)/pi+eps])
    g = animatedline; % subplot(2,1,2)


    for n2Time = 2:10:length(TimeVec)
        figure(f1)
        subplot(2,1,1)
        hold on
        addpoints(hprevoius,xBatPos(max(1,n2Time-20)),yBatPos(max(1,n2Time-20)));
        drawnow
        addpoints(h,xBatPos(n2Time),yBatPos(n2Time));
        drawnow 
        subplot(2,1,2)
        addpoints(g,TimeVec(n2Time),Teta(n2Time)/pi);
        drawnow
    end
end % AnimationFlag
end % BatFlight Fuction :-)
    
    





