function [IsAnyPhantomDetectedFlag, PhantomDetectionsStrct] = IsPhantomEchoesDetectected(InterEchoesStrct, BAT, RxBatNum, AllParams)

% This function calcultes whether the echoes from conspecifics signals are detected as targets
[ThMaxRx, IndEcho] = find([InterEchoesStrct.MaxRxPowerDB] >= ...
    AllParams.BatSonarParams.PulseDetectionTH);
if isempty(ThMaxRx)
    IsAnyPhantomDetectedFlag = 0;
    PhantomDetectionsStrct = [];
else % if isempty(ThMaxRx)
    PhantomDetectionsStrct = BAT(RxBatNum).PhantomEchoesStrct;
    IsAnyPhantomDetectedFlag =1;
    for kEcho = IndEcho
        kk=1;
        for kTimes = ...
                [InterEchoesStrct(kEcho).StartnTime : InterEchoesStrct(kEcho).EndnTime]
            
            PhantomDetectionsStrct(kTimes).nTime = [PhantomDetectionsStrct(kTimes).nTime, kTimes];
            PhantomDetectionsStrct(kTimes).PreyNum = [PhantomDetectionsStrct(kTimes).PreyNum , ...
                InterEchoesStrct(kEcho).PreyNum ] ;
            PhantomDetectionsStrct(kTimes).TxBatNum = [PhantomDetectionsStrct(kTimes).TxBatNum , ...
                InterEchoesStrct(kEcho).TxBatNum ] ;
            PhantomDetectionsStrct(kTimes).SignalPhase = [PhantomDetectionsStrct(kTimes).SignalPhase , ...
                {InterEchoesStrct(kEcho).SignalPhase} ] ;
            PhantomDetectionsStrct(kTimes).StartnTime = [PhantomDetectionsStrct(kTimes).StartnTime , ...
                InterEchoesStrct(kEcho).StartnTime ] ;
            PhantomDetectionsStrct(kTimes).EndnTime = [PhantomDetectionsStrct(kTimes).EndnTime , ...
                InterEchoesStrct(kEcho).EndnTime ] ;
            PhantomDetectionsStrct(kTimes).Freqs = [PhantomDetectionsStrct(kTimes).Freqs , ...
                InterEchoesStrct(kEcho).Freqs(kk) ] ;
            PhantomDetectionsStrct(kTimes).RxPowersDB = [PhantomDetectionsStrct(kTimes).RxPowersDB , ...
                InterEchoesStrct(kEcho).RxPowersDB(kk) ] ;
            PhantomDetectionsStrct(kTimes).MaxRxPowerDB = [PhantomDetectionsStrct(kTimes).MaxRxPowerDB , ...
                InterEchoesStrct(kEcho).MaxRxPowerDB ] ;
            PhantomDetectionsStrct(kTimes).Distance2PhantomReal = [PhantomDetectionsStrct(kTimes).Distance2PhantomReal , ...
                InterEchoesStrct(kEcho).Distance2PhantomReal ] ;
            PhantomDetectionsStrct(kTimes).Angle2Phantom = [PhantomDetectionsStrct(kTimes).Angle2Phantom , ...
                InterEchoesStrct(kEcho).Angle2Phantom ] ;
            kk = kk+1;
            
            
        end % for kTimes
    end % for kEcho
end %if isempty(ThMaxRx)