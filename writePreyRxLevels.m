function  [RxnTimes, RxLvls, RxBatNumTx, RxBatPulseNum]  = writePreyRxLevels(RxLvlStrct, ObjStruct, kPrey, kTx, PulseNum, BatOrPrey)

% this function write the rxlevels of the ecolocation calls reiveved by the
% prey or the rxlevels of the jamming calls received  by the bats
%
% inputs:
    % RxLvlStrct: the struct of rx levels and timings (output ofth efunction calcSigRxLvl)
    % ObjStruct: the struct of the receivers[ BatDATA.BAT(BatNum) or BatDATA.PREY(kPrey)]
    % kPrey - the index of the prey
    % kTx- the index of the transmitting (bat or k)
    % PulseNum - the number of the pulse
    % BatOrPrey - the rx object 'Bat' | 'Prey' (default is 'Prey')



%%%% 
% insert the recent data by order of the time of recption
if nargin < 6
    BatOrPrey = 'Prey';
end

switch BatOrPrey
    case 'Bat'
        RxnTimesF      = 'JammingPreyRxnTimes';
        RxLvlsF        = 'JammingPreyRxLvls';
        RxBatNumTxF    = 'JammingPreyNumTx';
        RxBatPulseNumF = 'JammingPreyPulseNum';
        curPrey_ix    = 1; % Only One Target and One Bat

    case 'Prey'
        RxnTimesF      = 'RxnTimes';
        RxLvlsF        = 'RxLvls';
        RxBatNumTxF    = 'RxBatNumTx';
        RxBatPulseNumF = 'RxBatPulseNum';
        curPrey_ix = find(RxLvlStrct.TargetIndex == kPrey);

end % switch

ntot = size(ObjStruct.(RxnTimesF),1);
ix = find(ObjStruct.(RxnTimesF) > RxLvlStrct.EchosTimes(curPrey_ix), 1);

if isempty(ix) 
    RxnTimes = [ObjStruct.(RxnTimesF), RxLvlStrct.EchosTimes(curPrey_ix)];
    RxLvls = [ObjStruct.(RxLvlsF), RxLvlStrct.RxLvl(curPrey_ix)];
    RxBatNumTx = [ObjStruct.(RxBatNumTxF), kTx];  
    RxBatPulseNum = [ObjStruct.(RxBatPulseNumF), PulseNum];
else  
    RxnTimes = [ObjStruct.(RxnTimesF)(1:ix-1), RxLvlStrct.EchosTimes(curPrey_ix), ObjStruct.(RxnTimesF)(1:ntot)];
    RxLvls = [ObjStruct.(RxLvlsF)(1:ix-1), RxLvlStrct.RxLvl(curPrey_ix), ObjStruct.(RxLvlsF)(1:ntot)];
    RxBatNumTx = [ObjStruct.(RxBatNumTxF)(1:ix-1), kTx, ObjStruct.(RxBatNumTxF)(1:ntot)];
    RxBatPulseNum = [ObjStruct.(RxBatPulseNumF)(1:ix-1), PulseNum, ObjStruct.(RxBatPulseNumF)(1:ntot)];
end % if iempty(ix) 
