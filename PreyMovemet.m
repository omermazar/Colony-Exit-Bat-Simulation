function [xPNext,yPNext, PTetaNext , PDtetaNext, CurrentVelocity, PObsAvoidManFlag] = ...
        PreyMovemet( xP, yP, TetaP, CurrentVelocity, CurrentDteta,...
        SimParameters, Terrain, ChangeDirCommand)

% [xPNext,yPNext,TetaPNext, ObsAvoidManFlag, DtetaNext, nTimeToMan] = ...
%         PreyMovemet(xP,yP, TetaP, CurrentVelocity, CurrentDteta,...
%         Obsticle, SimParameters)


% This funtion will calculate tne new position(x,y,teta) of the Pry
% output arguments: 
    % xPNext -the x coordinate ,  yPNext -the y coordinate , 
    % TetaPNext- the angle of the prey flight direction
    % ObsAvoidManFlag - wether The prey is manuvering to avoid Obs, should
    % be 0 (zero)
    % nTimeToMan - the period time needed to this manuver (normaly 1) 
% Inptus: 
    % xP,yP, TetaP- current position
    % CurrentVelocity = the velocity of the Prey
    % FlightType = 'linear', 'lircle', random'
    % CurrentVelocity- the CurrentVelocity of the prey
    % CurrentDteta =- the differntial angle of the prey ' from its accelartion
    % SimParameters - a Struct with the foolwing tstrcut:
    %       SimParameters.SimParamsB - xyResolution and SampleTime
    %       SimParameters.BatFlightParameters - MaxAccel, FlightType
    %       SimParameters.PreyFlightParameters     - ...
    % ChangeDirCommand - '1' - if need to cange, '0' to star with the
    %       current dteta
    % Modified by Omer 9/10/2017
    
    % Set Paramters from Input
    % Vladimir
    %%%%% Nov2021
    
    FlightType = SimParameters.PreyFlightParams.PFlightTypeFlag;
    xyResolution = SimParameters.SimParams.xyResolution;
    SampleTime = SimParameters.SimParams.SampleTime;
    MaxAccel = SimParameters.PreyFlightParams.PMaxAccelaration * SampleTime^2/xyResolution;
    NominalAccel = SimParameters.PreyFlightParams.PNominalAccel * SampleTime^2/xyResolution;
    NominalRadius = CurrentVelocity.^2/NominalAccel;
    MinRadius =  CurrentVelocity.^2/MaxAccel;
    PmaxVelocity = SimParameters.PreyFlightParams.PNominalVelocity * SampleTime/xyResolution * 1.2;
    PminVelocity = SimParameters.PreyFlightParams.PNominalVelocity *SampleTime/xyResolution * 0.8;
%     ObsDist = Obsticle.Distance;
%     ObsAngle = Obsticle.Alpha;
    StdRand = sqrt(2)*SimParameters.PreyFlightParams.PStdRandom ;
%   IsObs = Obsticle.IsObs;
    
    minDistAllowed = max(SimParameters.PreyFlightParams.PMinDistanceAllowed/xyResolution, CurrentVelocity/SampleTime * 0.15);
  
    
    
    KdTeta = 1;  % Parameter for circleflight
    MaxDTetaPrey = MaxAccel/CurrentVelocity; % Maximum Radian Velocity [rad/sec]
    NominalOmega = NominalAccel/CurrentVelocity;
%     MaxDTetaPrey = MaxOmega; % the maximum change in time-frame
    NominalDTetaPrey = NominalOmega;
    RegTimeToDecision = 1; % The Prey has no SOnar (it sees the obsticle regular- in free terrain
    minTimeToDecision = 1; % PRF when h
    
    % OutPuts INIT
    PreyManFlag = 0;
    PObsAvoidManFlag = 0;
    nTimeToMan = 0;
    
    kNearestObj = 1; % the room limits
    % if it is time to make the decision, make it
    % time to make decisin = sonar pulse transmitted (recived...?)
    
    
    
%     %% Make decision
%     
%     %Initialization
%     ObsStruct.Distances = [];
%     ObsStruct.TargetAngle = [];
%     sz = size(Terrain);
%     sz=sz(1);
%     IsObs = 0; 
%     
%     % Room Limits
%     % [in, on] = inpolygon(xP, yP, Terrain(1,:), Terrain(2,:));
%     %  lntersecId = find(in==0,1,'first');
%     d  = abs(p_poly_dist(xP, yP, Terrain(1,:), Terrain(2,:) ) ) ;
%     lntersecId = d <= minDistAllowed;
%     if lntersecId>0
%         IsObs = 1;
%     end
%     
%     % Intersection objects
%     if ~strcmp(SimParameters.TerrainParams.Type,  'Empty Room') && ... 
%         ~strcmp(SimParameters.TerrainParams.Type,  'Free Space')
%         for i=3:+2:sz-1 %for every object loop
%             % [in, on] = inpolygon(xP, yP, Terrain(i,:), Terrain(i+1,:));
%             % lntersecId = find(in,1,'first');
%             d  = abs(p_poly_dist(xP, yP, Terrain(i,:), Terrain(i+1,:) ) ) ;
%             lntersecId = d <= minDistAllowed;
%             if lntersecId>0
%                 kNearestObj = i;
%                 IsObs = 1;
%             end
%         end
%     end % if ~strcmp(SimParameters.TerrainParams.Type,  'Empty Room')
    
    %% Movment
    
    if ChangeDirCommand == 1
        
         %% Make decision
    
    %Initialization
    ObsStruct.Distances = [];
    ObsStruct.TargetAngle = [];
    sz = size(Terrain);
    sz=sz(1);
    IsObs = 0; 
    
    % Room Limits
    % [in, on] = inpolygon(xP, yP, Terrain(1,:), Terrain(2,:));
    %  lntersecId = find(in==0,1,'first');
    d  = abs(p_poly_dist(xP, yP, Terrain(1,:), Terrain(2,:) ) ) ;
    lntersecId = d <= minDistAllowed;
    if lntersecId>0
        IsObs = 1;
    end
    
    % Intersection objects
    if ~strcmp(SimParameters.TerrainParams.Type,  'Empty Room') && ... 
        ~strcmp(SimParameters.TerrainParams.Type,  'Free Space')
        for i=3:+2:sz-1 %for every object loop
            % [in, on] = inpolygon(xP, yP, Terrain(i,:), Terrain(i+1,:));
            % lntersecId = find(in,1,'first');
            d  = abs(p_poly_dist(xP, yP, Terrain(i,:), Terrain(i+1,:) ) ) ;
            lntersecId = d <= minDistAllowed;
            if lntersecId>0
                kNearestObj = i;
                IsObs = 1;
            end
        end
    end % if ~strcmp(SimParameters.TerrainParams.Type,  'Empty Room')

        if IsObs == 0 % no need to manuver away from obsticle
       
            switch FlightType
                case 'linear'
                    PDtetaNext = 0; 
                case 'circle'
                    PDtetaNext =  NominalDTetaPrey*KdTeta;
                case 'random'
                   CurrentVelocity = CurrentVelocity + 0.1*randn(1,1);
                   % keeping the velocity between PminVelocity and  PmaxVelocity
                   if CurrentVelocity >=  PmaxVelocity
                       CurrentVelocity = PmaxVelocity;
                   elseif CurrentVelocity <=  PminVelocity
                       CurrentVelocity = PminVelocity;
                   end % if CurrentVelocity > =  PmaxVelocity
                   PDtetaNext = MaxDTetaPrey*StdRand*randn;
                case 'static'
                     CurrentVelocity = 0;
                     PDtetaNext = 0;
                otherwise
                    disp('wrong Flight Type Input')
            end % switch FlightType
% % % %         elseif ChangeDirCommand == 0 % if ChangeDirCommand == 1
% % % %             PDtetaNext = CurrentDteta;
% % % %         end %if ChangeDirCommand == 1
        
        PTetaNext = TetaP + PDtetaNext;

        elseif  IsObs == 1

            % PDtetaNext = -2*TetaP; % pi/2-TetaP;
            % if PDtetaNext == 0
            %    PDtetaNext = pi;
            % end

            %%%%% Get away from the Wall %%%
            [~, x_poly, y_poly] = p_poly_dist(xP, yP, Terrain(kNearestObj,:), Terrain(kNearestObj+1,:) );
            PTetaNext  = atan2(yP - y_poly , xP - x_poly);
            PDtetaNext = 0; % PTetaNext - TetaP;

        end % if ObsFlad

    elseif ChangeDirCommand == 0 % if ChangeDirCommand == 1
        PDtetaNext = CurrentDteta;
        PTetaNext   = TetaP + PDtetaNext;
    end %if ChangeDirCommand == 1
          
    %[-pi,pi]
    % PTetaNext = pi*sawtooth(PTetaNext-pi);
    PTetaNext = PTetaNext -(2*pi)*floor((PTetaNext+pi)/(2*pi));
    
    xPNext = xP + CurrentVelocity.*cos(PTetaNext) ;
    yPNext = yP + CurrentVelocity.*sin(PTetaNext) ;

  end% PreyMovement Function