function [finds_struct] = update_Consps_FindsStruct(BAT)

% This function updsates the structure odf findings for conspeficifcs
% inputs: BAT = BatDATA.BAT(#)

% init
CurrPulseNum = BAT.CurrPulseNum;
% finds_struct = BAT.Consps_FindsStruct(1);
finds_struct = BAT.Consps_FindsStruct(CurrPulseNum); % empty 

% update
Detected_Consps = BAT.Consps_MaskingStruct(CurrPulseNum).DetectedPreys;
[~, ~, Detected_idx] = intersect(Detected_Consps, BAT.OtherBatsPolarCoordinates(CurrPulseNum).TargetsID);
[Detected_Consps_UnMasked , TargetInd] = ...
                                setdiff( Detected_Consps ,...
                                BAT.Consps_MaskingStruct(CurrPulseNum).MaskedPreys ,'stable'); 
% PreyIndex = BAT.DetectedPreysUnMasked;

finds_struct.TransmittedPulseNum = CurrPulseNum;
finds_struct.IsAnyPreyDetected = numel(Detected_Consps) > 0;
finds_struct.DetectedPreyNum = Detected_Consps_UnMasked;
finds_struct.DetectedTimes = BAT.Detected_Consps_Times(TargetInd);
finds_struct.PreyNumToHunt = 0;

finds_struct.DetecectedPreyWithOutInterference = Detected_Consps;
finds_struct.DetectedTimesWithOutInterfernece = BAT.Detected_Consps_Times;
finds_struct.MissedPreysByProbabilty = [];

finds_struct.Dist2DetectedPrey = ...
    BAT.OtherBatsPolarCoordinates(CurrPulseNum).Distances(Detected_idx);
finds_struct.Angle2DetectedPrey = ...
    BAT.OtherBatsPolarCoordinates(CurrPulseNum).TargetAngle(Detected_idx);
% finds_struct.PreyRelativeDirection = ...
%     BAT.OtherBatsPolarCoordinates(CurrPulseNum).Bat2TargetRelativeAngle(Detected_idx);

if finds_struct.IsAnyPreyDetected
    finds_struct.MaskedPreys = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).MaskedPreys;
    finds_struct.RxPowerOfDetectedPreys = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).DetectedPreysRxPower;
    finds_struct.SIROfDetectedPreys = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).DetecetedPrey2InterferenceRatioDB;
    
    finds_struct.SelfCorrealationMaxdB = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).SelfCorrealationMaxdB;
    finds_struct.InterCorrelationMaxdB = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).InterCorrelationMaxdB;
    
    finds_struct.xFinds = [];
    finds_struct.yFinds = [];
    
    finds_struct.DetectionsDFerror = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).DetectionsDFerror;
    finds_struct.DetectionsRangeErr = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).DetectionsRangeErr;
    finds_struct.DetectionsRelativeDirectionErr = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).DetectionsRelativeDirectionErr;
    
    % FilterBank
    finds_struct.FB_unmasked_prey = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).FB_unmasked_prey;
    finds_struct.FB_unmasked_delays = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).FB_unmasked_delays;
    finds_struct.FB_detected_masking_delays = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).FB_detected_masking_delays;
    finds_struct.FB_estimated_masker = ...
        BAT.Consps_MaskingStruct(CurrPulseNum).FB_estimated_masker;
    
    % Estmateed location: (New Sep2022) 
    [~,ixDetConsps] = ismember(finds_struct.DetectedPreyNum, ...
                                finds_struct.DetecectedPreyWithOutInterference );              
    finds_struct.Distances = finds_struct.Dist2DetectedPrey(ixDetConsps)  + finds_struct.DetectionsRangeErr(ixDetConsps);
    finds_struct.Angles    = finds_struct.Angle2DetectedPrey(ixDetConsps) + finds_struct.DetectionsDFerror(ixDetConsps);

    
else % Update nonhunting
   finds_struct.PreyNumToHunt = ...
        0;
    finds_struct.IsHuntedPreyMasked = ...
        0;
end %if finds_struct.IsAnyPreyDetected

% Clutter
finds_struct.ClutteredPrey = ...
    BAT.ClutteredPreyVec;
finds_struct.Clutter_nTimes = ...
    BAT.ClutterednTimesVec;

