function [finds_struct] = updateFindsStruct(BAT, CurrPulseNum, nTime, detectedObj, varargin)

% This function updsates the structure odf findings for conspeficifcs
% inputs:   BAT = BatDATA.BAT(#)
%           detectedObj = 'Prey' | 'Obs' | 'Consps'

%%%% New June2024 - Confusion between own echoes and ehcoes from conspecifics calls
% varagin: 'confusionFlag' - if the bat can discriminate baetween own echoes
%                            and echoes from cons, the value is false (default)
%          'echoTable'     - the table of the relavant echoes ecived bt the
%                            bat. echoTable = BAT.echoesFullTable. NOTICE:
%                            the table should be filtered outside for the
%                            relevant echo-types.
% if confusionFlag is true - the referncens targetts are taken from
% BAT.echoesFullTable insteats of pLoc_struct    = BAT.ObsInBeamStruct(CurrPulseNum);

%% Input- Changed Jam Moth Apr2022
if nargin < 12
    debug_mode = 0;
end

if nargin < 13
    DetectedObj = 'Prey';
end

% default params 
confusionFlag = false; % default - bat can discriminate baetween own echoes and echoes from cons
echoTable = table();

if ~isempty(varargin)
    if mod(numel(varargin), 2) ~= 0
        warning('wrong input varargin')
    else
        % inParam = cell(numel(varargin),1);
        for k = 1:2:numel(varargin)
            % option 1
            % inParam.varagin{k} = varagin{k+1};
            % option 2
            switch varargin{k}
                case 'confusionFlag'
                    confusionFlag = varargin{k+1};
                        
                case 'echoTable'
                    echoTable = varargin{k+1};
                
                case 'AllParams'
                    AllParams = varargin{k+1};;

                otherwise
                    disp('unkwnown input parameter')
            end % switch
        end % for k
    end % if mod
end % if ~isempty(varargin)

% init
% CurrPulseNum = BAT.CurrPulseNum;
% finds_struct = BAT.Consps_FindsStruct(1);

switch detectedObj
    case 'Prey' 
        finds_struct   = BAT.PreyFindsStruct(CurrPulseNum); % empty 
        masking_struct = BAT.FindsMaskingStruct(CurrPulseNum);
        pLoc_struct    = BAT.CurrDetectedPreysPolarCoordinate;
        [~,inBeamIx]   = ismember(masking_struct.DetectedPreys, BAT.EchosFromPreyStruct(CurrPulseNum).TargetIndex);
        rxTimes        = BAT.EchosFromPreyStruct(CurrPulseNum).EchosTimes(inBeamIx) + BAT.TransmittedPulsesStruct(CurrPulseNum).StartPulseTime;
    
    case 'Obs' 
        finds_struct   = BAT.ObsFindsStruct(CurrPulseNum+1); % empty
        masking_struct = BAT.Obs_MaskingStruct(CurrPulseNum);
        pLoc_struct    = BAT.ObsInBeamStruct(CurrPulseNum);
        
        pLoc_struct.NumOfTargets = numel(pLoc_struct.DetectetedTargetsVec);
        pLoc_struct.TargetsID    = pLoc_struct.DetectetedTargetsVec;
        if ~confusionFlag % New June2024
            [~,inBeamIx]   = ismember(masking_struct.DetectedPreys, BAT.EchosFromObsStruct(CurrPulseNum).TargetIndex);
            inBeamIx       = nonzeros(inBeamIx); % New June2024 - 
            rxTimes        = BAT.EchosFromObsStruct(CurrPulseNum).EchosTimes(inBeamIx) + BAT.TransmittedPulsesStruct(CurrPulseNum).StartPulseTime;
        end % if ~confusionFlag

    case 'Consps'
        finds_struct   = BAT.Consps_FindsStruct(CurrPulseNum); % empty 
        masking_struct = BAT.Consps_MaskingStruct(CurrPulseNum); 
        pLoc_struct    = BAT.OtherBatsPolarCoordinates(CurrPulseNum);
        [~,inBeamIx]   = ismember(masking_struct.DetectedPreys, BAT.EchosFromConspStruct(CurrPulseNum).TargetIndex);
        rxTimes        = BAT.EchosFromConspStruct(CurrPulseNum).EchosTimes(inBeamIx) + BAT.TransmittedPulsesStruct(CurrPulseNum).StartPulseTime;
        
end % switch

% The position of the Tx bat
calcTime = BAT.TransmittedPulsesStruct(CurrPulseNum).StartPulseTime;
xBat           = BAT.xBati(calcTime);
yBat           = BAT.yBati(calcTime);
tetaBat        = BAT.Teta(calcTime);

%%%% if confuisonFlag tha override pLoc_sttruct New June2024 
if confusionFlag

    SampleTime = AllParams.SimParams.SampleTime;

    TxTime = BAT.TransmittedPulsesStruct(CurrPulseNum).StartPulseTime*SampleTime;
    maxTime = TxTime+ BAT.TransmittedPulsesStruct(CurrPulseNum).IPItoNextPulse*SampleTime;
    % the relevant signals are obstacle in the cuurent time of pulse:
    ix = (BAT.echoesFullTable.time >= TxTime) & (BAT.echoesFullTable.time <= maxTime) & ...
        (BAT.echoesFullTable.echoType == "ownObs" | BAT.echoesFullTable.echoType == "ConspsObsEchoes");
    echoTable = BAT.echoesFullTable(ix,:);
    
    % recontruct the pLoc
    %%% The estimated distancesare taken from the unexpected time of
    %%% recption 
    pLoc_struct.Distances   = (echoTable.time' - TxTime)*(343/2) / AllParams.SimParams.xyResolution;
    % % % % pLoc_struct.Distances   = echoTable.distance'/ AllParams.SimParams.xyResolution;
    pLoc_struct.TargetAngle = echoTable.directionOfArrival';
    pLoc_struct.xOBSTC      = xBat + pLoc_struct.Distances .* cos(pLoc_struct.TargetAngle+tetaBat);
    pLoc_struct.yOBSTC      = yBat + pLoc_struct.Distances .* sin(pLoc_struct.TargetAngle+tetaBat);
    pLoc_struct.DetectetedTargetsVec = echoTable.targetIndex';
    pLoc_struct.NumOfTargets = numel(pLoc_struct.DetectetedTargetsVec);
    pLoc_struct.TargetsID    = pLoc_struct.DetectetedTargetsVec;

    % [~,inBeamIx]   = ismember(masking_struct.DetectedPreys, echoTable.targetIndex);
    % inBeamIx       = nonzeros(inBeamIx); % New June2024 -
   
    rxTimes        = echoTable.time' / SampleTime;
    targetInd_pLocToMaskStruct = find(ismember(pLoc_struct.TargetsID, masking_struct.DetectedPreys));
end


% update
detectedTargets = masking_struct.DetectedPreys;
[~, ~, Detected_idx] = intersect(detectedTargets, pLoc_struct.TargetsID,'stable');
[detectedTargets_UnMasked , TargetInd] = setdiff( detectedTargets, masking_struct.MaskedPreys ,'stable');
ixUnMasked_in_ploc_struct = Detected_idx(TargetInd);
% PreyIndex = BAT.DetectedPreysUnMasked;

finds_struct.TransmittedPulseNum = CurrPulseNum;
finds_struct.IsAnyPreyDetected   = numel(detectedTargets) > 0;
finds_struct.DetectedPreyNum     = detectedTargets_UnMasked;
finds_struct.PreyNumToHunt       = 0;

finds_struct.DetecectedPreyWithOutInterference = detectedTargets;
if ~confusionFlag
    finds_struct.DetectedTimes       = rxTimes(TargetInd);
    finds_struct.DetectedTimesWithOutInterfernece  = rxTimes;
else
    finds_struct.DetectedTimes       = rxTimes(ixUnMasked_in_ploc_struct);
    finds_struct.DetectedTimesWithOutInterfernece  = rxTimes(Detected_idx);
end
finds_struct.MissedPreysByProbabilty = [];

finds_struct.Dist2DetectedPrey  = pLoc_struct.Distances(ixUnMasked_in_ploc_struct); % Detected_idx
finds_struct.Angle2DetectedPrey = pLoc_struct.TargetAngle(ixUnMasked_in_ploc_struct); % Detected_idx

if finds_struct.IsAnyPreyDetected
    finds_struct.MaskedPreys            = masking_struct.MaskedPreys;
    finds_struct.RxPowerOfDetectedPreys = masking_struct.DetectedPreysRxPower; % Cahnge Aug2024
    finds_struct.SIROfDetectedPreys     = masking_struct.DetecetedPrey2InterferenceRatioDB(TargetInd);
    finds_struct.SelfCorrealationMaxdB  = masking_struct.SelfCorrealationMaxdB(TargetInd);
    finds_struct.InterCorrelationMaxdB  = masking_struct.InterCorrelationMaxdB(TargetInd);
    
    finds_struct.DetectionsDFerror  =             masking_struct.DetectionsDFerror(TargetInd);
    finds_struct.DetectionsRangeErr =             masking_struct.DetectionsRangeErr(TargetInd);
    finds_struct.DetectionsRelativeDirectionErr = masking_struct.DetectionsRelativeDirectionErr(TargetInd);

    switch detectedObj 
        case 'Obs'
            finds_struct.xFinds = pLoc_struct.xOBSTC(ixUnMasked_in_ploc_struct);
            finds_struct.yFinds = pLoc_struct.yOBSTC(ixUnMasked_in_ploc_struct);
        case 'Prey'
            finds_struct.Bat2TargetRelativeAngle = pLoc_struct.Bat2TargetRelativeAngle(ixUnMasked_in_ploc_struct) + ...
                finds_struct.DetectionsRelativeDirectionErr;
            finds_struct.xFinds = BAT.xFindsCurrentVec(ixUnMasked_in_ploc_struct);
            finds_struct.yFinds = BAT.yFindsCurrentVec(ixUnMasked_in_ploc_struct);
        otherwise %if strcmp(detectedObj, 'Obs')
            finds_struct.Bat2TargetRelativeAngle = pLoc_struct.Bat2TargetRelativeAngle(ixUnMasked_in_ploc_struct) + ...
                finds_struct.DetectionsRelativeDirectionErr;
            finds_struct.xFinds = xBat + finds_struct.Dist2DetectedPrey.*cos(tetaBat + finds_struct.Angle2DetectedPrey);
            finds_struct.yFinds = yBat + finds_struct.Dist2DetectedPrey.*sin(tetaBat + finds_struct.Angle2DetectedPrey);

    end % switch detectedObj 
    
    % Estimated Localization
    finds_struct.Distances = finds_struct.Dist2DetectedPrey  + finds_struct.DetectionsRangeErr;
    finds_struct.Angles    = finds_struct.Angle2DetectedPrey + finds_struct.DetectionsDFerror;

    % The estimated [x,y] position
    finds_struct.xFindsEst = xBat + finds_struct.Distances.*cos(finds_struct.Angles+tetaBat);
    finds_struct.yFindsEst = yBat + finds_struct.Distances.*sin(finds_struct.Angles+tetaBat);

else % Update nonhunting
   finds_struct.PreyNumToHunt       = 0;
   finds_struct.IsHuntedPreyMasked =  0;
end %if finds_struct.IsAnyPreyDetected

% Clutter
finds_struct.ClutteredPrey =  BAT.ClutteredPreyVec;
finds_struct.Clutter_nTimes = BAT.ClutterednTimesVec;

