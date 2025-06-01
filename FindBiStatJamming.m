function [IsJammed] = ...
    FindBiStatJamming(CurrentBiStatPulses, InterferenceFullStruct, AllParams)

% this function calculates whether the Bi-Stat detection is jammed
% the function checks the masking windooe in the certain frequency band of the bi-stat detection\\

 % from mili-sec to samples
MaskingFwdnTime = AllParams.BatSonarParams.MaskingFwdTime / AllParams.SimParams.SampleTime *1e-3 ;
MaskingBckdnTime = AllParams.BatSonarParams.MaskingBckdTime / AllParams.SimParams.SampleTime *1e-3 ;
% the frerquency Band
HalfBandWidth = 0.5; % kHz
JamBand = [CurrentBiStatPulses.PreyToHuntFreq -  HalfBandWidth , ...
    CurrentBiStatPulses.PreyToHuntFreq +  HalfBandWidth];
StartnTime = max(1, CurrentBiStatPulses.PreyToHuntTime - MaskingFwdnTime+1 );
EndnTime = min( AllParams.SimParams.SimulationTime / AllParams.SimParams.SampleTime , ...
    CurrentBiStatPulses.PreyToHuntTime + MaskingBckdnTime);

IntFreqs = [InterferenceFullStruct(StartnTime:EndnTime).Freqs];
IntPowers = [InterferenceFullStruct(StartnTime:EndnTime).Power];
% the times the innterference freqs are in the band
MaskFreqIdx =   IntFreqs >= JamBand(1) & IntFreqs <= JamBand(2);
MaskingPowersDB = 10*log10(abs(IntPowers(MaskFreqIdx)));

IsJammed = sum(MaskingPowersDB > ...
    (CurrentBiStatPulses.PreyToHuntRxPowerDB + AllParams.BatSonarParams.MaskingBckdSIRth) ) > 0;




                            