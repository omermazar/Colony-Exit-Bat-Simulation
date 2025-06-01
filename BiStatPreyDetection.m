function [CurrPulseBiStat ] = ...
    BiStatPreyDetection(BiStatDetection, TransmittedPulsesStruct, PreyFindsStruct, AllParams)
% check  only for undetected prey if the bi-stat echo is stronger than detecton TH 
% and than check if it would change the bat decision - only if it's stronger than the srorngest echo from prey 

NumOfPreyItems = AllParams.SimParams.TotalPreysNumber;
NoiseLeveldB = AllParams.BatSonarParams.NoiseLeveldB;
PulseDetectionTH = AllParams.BatSonarParams.PulseDetectionTH;

PreyToHuntBiSonar = 0;
PreyToHuntRxPowerDB = 0;
PreyToHuntTime = 0;
PreyTohHuntFreq =0;
Dist2HuntedPrey = -1;
Angle2HuntedPrey = pi;

minTime = round(TransmittedPulsesStruct.StartPulseTime + TransmittedPulsesStruct.PulseWidth);
maxTime = TransmittedPulsesStruct.StartPulseTime +   ...
    max( TransmittedPulsesStruct.IPItoNextPulse , TransmittedPulsesStruct.PulseWidth+1 );
ReceptionTimes = minTime  : maxTime; 

BiSonarCurrentRxPowers = BiStatDetection.AllPreyRxPowerMat(:,ReceptionTimes);
[MaxPowers, Delays] = max( BiSonarCurrentRxPowers,[],2);
UnDetectedPrey = setdiff(1:NumOfPreyItems, PreyFindsStruct.DetectedPreyNum);
MaxPowersDB = 10*log10(abs(MaxPowers));
PreyDetectedBiSonar = find(MaxPowersDB > NoiseLeveldB + PulseDetectionTH);
% DetectedBiSonarPowerDB = MaxPowersDB(PreyDetectedBiSonar);
BiSonarTimes = minTime + Delays- 1;
PreyDetectedOnlyByBiSonar = intersect(UnDetectedPrey, PreyDetectedBiSonar);

[MaxRx, idx] = max(MaxPowers(PreyDetectedOnlyByBiSonar));

if MaxRx > 10^(max(PreyFindsStruct.RxPowerOfDetectedPreys)/10)
    PreyToHuntBiSonar = PreyDetectedOnlyByBiSonar(idx);
    PreyToHuntRxPowerDB = 10*log10(abs(MaxRx));
    PreyToHuntTime = BiSonarTimes(PreyToHuntBiSonar);
    PreyTohHuntFreq = BiStatDetection.AllPFreqMat(PreyToHuntBiSonar, PreyToHuntTime);
    Dist2HuntedPrey = BiStatDetection.Prey(PreyToHuntBiSonar).Dist2Prey(PreyToHuntTime);
    Angle2HuntedPrey = BiStatDetection.Prey(PreyToHuntBiSonar).Angle2Prey(PreyToHuntTime);
%     BiSonarCurrentRxFreqs = BiStatDetection.AllPFreqMat(:,ReceptionTimes);
%     PreyTohHuntFreq = BiSonarCurrentRxFreqs(PreyToHuntBiSonar, PreyToHuntTime);
end % if MaxRx > max(RxPowerOfDetectedPreys)

CurrPulseBiStat = struct(...
    'PreyDetectedBiSonar',  PreyDetectedBiSonar', ...
    'PreyDetectedOnlyByBiSonar',  PreyDetectedOnlyByBiSonar', ...
    'PreyToHuntBiSonar' , PreyToHuntBiSonar, ...
    'PreyToHuntRxPowerDB' , PreyToHuntRxPowerDB, ...
    'PreyToHuntTime', PreyToHuntTime, ...
    'PreyToHuntFreq', PreyTohHuntFreq, ...
    'Dist2HuntedPrey', Dist2HuntedPrey, ...
    'Angle2HuntedPrey', Angle2HuntedPrey, ...
    'IsPreyJammed', 0);
