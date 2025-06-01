function[Params]= CreateSimDefaultParams(varargin) 

% CtreateSimDefaultParams(FileOperation, PathName)
% This function creates, saves or updates the default parmetes for the bat flight simulation
% the functioמ can save the parametes into mat-file
% Inputs :
%           None : the function will return the strcucrure Params
%           FileOperation = 'save'  - to file
%                           'update' - update parmas and save
%           PathName = The Path needede to save
%           Updated params = the struct to update
% OutPuts :
%           Params : Strucutre ogf the parmeters needed for the simulation
%           Params.SimParams : the General parameters
%           Params.BatFlight : paramters for the flight of the Bat
%           Params.Sonar :  paramters of the sonar
% Example: CtreateSimDefaultParams('save','\\c:\MATLAB\...')

% Checking the INPUTS:
SaveFlag = 0;
if nargin ~= 0
    SaveFlag = 1;
    FileOperation = varargin{1};
    switch FileOperation
        case 'save'   
            switch nargin
                case 1
                  SavePath = 'C:\Users\DELL\Documents\עומר\לימודים\מחקר\MATLAB Bat Simulation\DATA\DefaultParameters';
                case 2
                   SavePath = varargin{2};
                otherwise
                    error('TOO MANY INPUTS'); 
            end % switch nargin
        case 'update'
            switch nargin
                case 2
                  SavePath = 'C:\Users\DELL\Documents\עומר\לימודים\מחקר\MATLAB Bat Simulation\DATA\DefaultParameters';
                  UpdatedParams = varargin{2};
                case 3
                   SavePath = varargin{2};
                   UpdatedParams = varargin{3};
                otherwise
                    error('TOO MANY INPUTS'); 
            end % switch nargin
        error('WRONG INPUTS');
    end % switch FileOperation
end % nargin ~= 0



%%%%%%%%
% Defaualt parameters for The Simulation
%%%%%%%%
% Simulation Paremeters
Params.SimParams.xyResolution = 0.01; % 0.01=1cm
Params.SimParams.SimulationTime = 20; % Sec
Params.SimParams.SampleTime = 1e-4; % Sec ,0.001= 1miliSec

Params.SimParams.AnimationFlag = true; % Animate Bat Movement in time
Params.SimParams.DebugFlag = 1; % run more DATA for inside functin debug
Params.SimParams.IsPreyFlag = 1; % 1 if there is a Prey
Params.SimParams.SoundV0 = 343;% the velocity of sound [m/sec]
Params.SimParams.TotalBatsNumber = 1; % the total number of bats
Params.SimParams.TotalPreysNumber = 1; % the total number of preys
%%%%%%%%%%%%%%%%%%%
% BAT FLIGHT - Basic Parameters
% nominal Velocity =1.5m/s
Params.BatFlight.MaxVelocity = 6 ; % m/s
Params.BatFlight.NominalVelocity = 1.5; % m/s 
Params.BatFlight.MaxAccelaration = 10; % [m/sec^2] ~2g 
Params.BatFlight.NominalAccel = 2; % m/s^2
Params.BatFlight.MinDistanceAllowed = 0.1; % [m] The minimum range the bat can be near Obsticle
Params.BatFlight.FlightTypeFlag = 'random'; % 'random' or 'linear' or 'circle'
Params.BatFlight.StdRandom = 1;   % The STD for Random Flight 
Params.BatFlight.PreyHuntingFlag = 1; % 1 if the Bat is hunting 
Params.BatFlight.minCatchPreyDistance = 0.05; %m THe Sistance the bat can vatch the prey

%%%%%%%%%%%%%%%%
% BAT SONAR BASIC
% PRF- max 100Hz, Nominal 20Hz
% Pulse_Time = max 3 , min 1 ms
% Beam = +-30deg
Params.BatSonar.BatBeamWidth = pi/6; % radians
%  the range the bat can detect obsicle
Params.BatSonar.BatDetectionRange = 3;% [m]  
Params.BatSonar.PulsePower = 1; % The Power of the Transmitted Pulse
%  the range that the pulse will interfere other bats
Params.BatSonar.BatInterferenceRange = Params.BatSonar.BatDetectionRange.^4;
Params.BatSonar.BatMaxPRF = Params.SimParams.SoundV0/ (2*Params.BatSonar.BatDetectionRange); %  Max PRF to Avoid distance ambiguity [1/sec] Trasmits per seconds
Params.BatSonar.BatNominalPRF = Params.BatSonar.BatMaxPRF/5; % [1/sec] Trasmits per seconds
Params.BatSonar.PulseTimeLong = 3*1e-3 ; % [sec] 
Params.BatSonar.PulseTimeShort = 1*1e-3; % [sec]
Params.BatSonar.BandWidth = 1; % [kHz] 
Params.BatSonar.CenterFreq = 80; % [1kHz]
Params.BatSonar.ChirpFlag = 1; % a flag for chirps 
Params.BatSonar.ChirpSpan = 4; % max freq - min frew of the chirp [kHz]
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% The Terrain Parameters- closed room 10X10X10m with Rectanguler Obsticle
% if Terrain(x,y)= 0 - there is no obsticle at x,y
% if Terrain(x,y)= 1 - there is an obsticle at x,y
% Room size: 10x10x10
Params.Terrain.Type = 'Empty Room';
Params.Terrain.Xmax = 10;
Params.Terrain.Ymax = 10; 
Params.Terrain.Zmax = 10; % meters
% 2 Rectangle Obsticles
Params.Terrain.OBbsXWidth = [0.5 , 2.5];
Params.Terrain.OBbsYWidth = [1 , 3]; % meters
Params.Terrain.OBbsDownLeftY = [2.5 , 5.5];
Params.Terrain.OBbsDownLeftX = [3, 6] ;% position of left buttom corner (meters)
% % Rectangle Obsticle #2
% Params.Terrain.OBbs2XWidth = 2.5;
% Params.Terrain.OBbs2YWidth = 3; % meters
% Params.Terrain.OBbs2DownLeftY = 5.5;
% Params.Terrain.OBbs2DownLeftX = 6;% position of left buttom corner (meters)

% PREY Flight Parameters
Params.PreyFlight.PMaxAccelaration = 2;
Params.PreyFlight.PmaxVelocity = 1; % m/s 
Params.PreyFlight.PNominalVelocity = 0.5; % m/s
Params.PreyFlight.PNominalAccel = 1; % m/s^2
Params.PreyFlight.PMinDistanceAllowed = 0.05; % [m] The minimum range the bat can be near Obsticle
Params.PreyFlight.PFlightTypeFlag = 'random'; % 'random' or 'linear' or 'circle'    
Params.PreyFlight.PStdRandom = 1;   % The STD for Random Flight 
Params.BatFlightParams.AvoidObsFlag = true; % to manuver away from absticles
% Outputs
if  SaveFlag == 1
    save(SavePath,'Params');
end
