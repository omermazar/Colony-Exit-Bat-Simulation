function [FindsAll] = obstacleDetectionSummary(BATx, AllParams, PreyObsFlag)

sampleTime   = AllParams.SimParams.SampleTime;
xyResolution = AllParams.SimParams.xyResolution;
nCalls       = BATx.NumOfTimesSonarTransmits;
switch PreyObsFlag
    case 'Prey'
        reqFindsStruct = BATx.PreyFindsStruct;
        MaskingStruct  = BATx.FindsMaskingStruct;
    case 'Obs'
        reqFindsStruct = BATx.ObsFindsStruct;
        MaskingStruct  = BATx.Obs_MaskingStruct;
    case 'Cosnps'
        reqFindsStruct = BATx.Consps_FindsStruct;
        MaskingStruct  = BATx.Consps_MaskingStruct;
end % switch

FindsAll(nCalls) = struct();
for k = 1:nCalls
    nTime = BATx.TransmittedPulsesStruct(k).StartPulseTime;
    FindsAll(k).PulseNum         = k;
    FindsAll(k).nTime            = nTime;
    FindsAll(k).t                = FindsAll(k).nTime*sampleTime;
    FindsAll(k).detectetObs      = reqFindsStruct(k).DetectedPreyNum;   
    FindsAll(k).maskedObs        = reqFindsStruct(k).MaskedPreys;  
    FindsAll(k).xObsReal         = reqFindsStruct(k).xFinds * xyResolution;  
    FindsAll(k).yObsReal         = reqFindsStruct(k).yFinds * xyResolution; 
    FindsAll(k).xObsEstimated    = (BATx.xBati(nTime) + reqFindsStruct(k).Distances .* cos(BATx.Teta(nTime) + reqFindsStruct(k).Angles)) * xyResolution;  
    FindsAll(k).yObsEstimated    = (BATx.yBati(nTime) + reqFindsStruct(k).Distances .* sin(BATx.Teta(nTime) + reqFindsStruct(k).Angles)) * xyResolution;
    FindsAll(k).Distances        = reqFindsStruct(k).Dist2DetectedPrey; 
    FindsAll(k).DistanceErrs     = reqFindsStruct(k).DetectionsRangeErr; 
    FindsAll(k).Angles           = reqFindsStruct(k).Angle2DetectedPrey; 
    FindsAll(k).AnglesErrs       = reqFindsStruct(k).DetectionsDFerror;
    FindsAll(k).RxLevel          = reqFindsStruct(k).RxPowerOfDetectedPreys;;
    FindsAll(k).SIR              = MaskingStruct(k).DetecetedPrey2InterferenceRatioDB;
   
end