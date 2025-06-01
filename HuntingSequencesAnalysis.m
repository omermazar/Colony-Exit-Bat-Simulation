function [Sequences] = HuntingSequencesAnalysis(BAT, AllParams)

% This guction unitedet the manuvers of hunting to long sequnces and
% summarizes the main outputs of each seqeunce

%%
%%% setup parameteres
AllHuntePreyVec = [BAT.ManueverCmdStruct.PreyNumToHunt]; 
DiffHundtedPreyVec= diff([0, AllHuntePreyVec]);
SampleTime = AllParams.SimParams.SampleTime;
xyResolution = AllParams.SimParams.xyResolution;
minCatchPreyDistance = AllParams.BatFlightParams.minCatchPreyDistance;
MinSequenceLength = 3; % the minimum length to star analyzing

% find the sequnces 
HuntManuversStartTimes = find([DiffHundtedPreyVec]);
HuntManueversLength = diff([HuntManuversStartTimes]);

IdxOfLongHunting = ( HuntManueversLength >= MinSequenceLength );
SequenceStarts = HuntManuversStartTimes(IdxOfLongHunting);
SequenceLengths = HuntManueversLength(IdxOfLongHunting) ;

HuntedPreyInSequnces = AllHuntePreyVec(SequenceStarts); % the Prey which the bat is hunting after
% Find non zeros Preys
IdNotzeros = (HuntedPreyInSequnces ~= 0);
HuntedPreyInSequnces = HuntedPreyInSequnces(IdNotzeros);
SequenceStarts = SequenceStarts(IdNotzeros);
SequenceLengths = SequenceLengths(IdNotzeros);
SequenceEnds = SequenceStarts + SequenceLengths-1;


% Unite sqeucences 
% if there are repeated sequences hunting the saame orey
% - than unit tem
MaxDiffFBetweenSequnces =5; % the maximpum number of pulses bewteen 2 repeated sequences
%  Repeated Sequnces
RepeatedSeqIdx = find(~diff([0 HuntedPreyInSequnces])); % index of o repeated seqeunces

% temp vars
UStarts= SequenceStarts;
UEnds = SequenceEnds;
if ~isempty(RepeatedSeqIdx)
    for k = RepeatedSeqIdx
        End2Start = SequenceStarts(k)-SequenceEnds(k-1);
        if End2Start <= MaxDiffFBetweenSequnces
            UStarts(k) =  UStarts(k-1);
            UEnds(k-1) = UEnds(k);
        end % if End2Start
    end % for k
    [UniteSequencesStarts, idx] = unique(UStarts);
    UniteSequenceEnds = UEnds(idx);
    UniteSequenceLengths = UniteSequenceEnds - UniteSequencesStarts +1;
    PreyInUniteSeqs = HuntedPreyInSequnces(idx);
    
else % if ~isempty(RepeatedSeqIdx)
    UniteSequencesStarts = SequenceStarts;
    UniteSequenceEnds = SequenceEnds;
    UniteSequenceLengths = SequenceLengths;
    PreyInUniteSeqs = HuntedPreyInSequnces;
end % if ~isempty(RepeatedSeqIdx)

%%% Filter out short sequences
MinSequenceLength = 4;

IdxOfLongSequnces = ( UniteSequenceLengths >= MinSequenceLength );
Sequences.NumOfHuntingSequnces = length(nonzeros(IdxOfLongSequnces));
Sequences.StartPulse = UniteSequencesStarts(IdxOfLongSequnces);
Sequences.LastPulse = UniteSequenceEnds(IdxOfLongSequnces);
Sequences.PreyInSeq = PreyInUniteSeqs(IdxOfLongSequnces);
Sequences.NumOfPulses = UniteSequenceLengths(IdxOfLongSequnces);
Sequences.StartTime = [BAT.TransmittedPulsesStruct([Sequences.StartPulse]).StartPulseTime]*SampleTime;
Sequences.EndTime = [BAT.TransmittedPulsesStruct([Sequences.LastPulse]).StartPulseTime]*SampleTime;
Sequences.Duration = [Sequences.EndTime] - [Sequences.StartTime];
Sequences.StartingStage = {BAT.ManueverCmdStruct([Sequences.StartPulse]).ManueverStage};
Sequences.EndStage =  {BAT.ManueverCmdStruct([Sequences.LastPulse-1]).ManueverStage};
Sequences.StartDistanceFromPrey = [BAT.ManueverCmdStruct([Sequences.StartPulse]).Dist2HuntedPrey]*xyResolution;
% dealing with  catches
MinEndDist = min([BAT.ManueverCmdStruct([Sequences.LastPulse-1]).Dist2HuntedPrey] ,...
    [BAT.ManueverCmdStruct([Sequences.LastPulse]).Dist2HuntedPrey])*xyResolution;
% MinEndDist = min([BAT.PreyFindsStruct([Sequences.LastPulse]).Dist2DetectedPrey])
Sequences.EndDistanceFromPrey = MinEndDist;
Sequences.IsCatch = MinEndDist <= 1.05*minCatchPreyDistance;

%  Masking to Sequence
Sequences.NumOfMakingToHunted = zeros(1, Sequences.NumOfHuntingSequnces );
% Sequences.EndDistanceFromPrey =  zeros(1, Sequences.NumOfHuntingSequnces );
for kk = 1:Sequences.NumOfHuntingSequnces 
    kkStartEnd = Sequences.StartPulse(kk):Sequences.LastPulse(kk);
    Sequences.NumOfMakingToHunted(kk) = sum(...
        [BAT.PreyFindsStruct(kkStartEnd).MaskedPreys] ==  Sequences.PreyInSeq(kk) ) ;
    
%     Hind = find(BAT.PreyFindsStruct(Sequences.LastPulse(kk)).DetecectedPreyWithOutInterference ...
%         == Sequences.PreyInSeq(kk) );
%     Sequences.EndDistanceFromPrey(kk) = ...
%         BAT.PreyFindsStruct(Sequences.LastPulse(kk)).Dist2DetectedPrey(Hind);
end % for kk
Sequences.MaskingRatio = Sequences.NumOfMakingToHunted ./ Sequences.NumOfPulses;
Sequences.NumOfCatches = sum(Sequences.IsCatch);
%%% OutPuts

% %%% plots
% figure
% hold on
% plot(AllHuntePreyVec);
% plot(HuntManuversStartTimes, AllHuntePreyVec(HuntManuversStartTimes),'o')
% plot(SequenceStarts, HuntedPreyInSequnces,'+k')
% plot(SequenceEnds, HuntedPreyInSequnces ,'+b')
% plot(UniteSequencesStarts, PreyInUniteSeqs, 'sm', 'MarkerSize', 8)
% plot(UniteSequenceEnds, PreyInUniteSeqs, 'sm', 'MarkerSize', 8)
