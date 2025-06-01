function [ IsPhantomDetected, PhantomAnalysisStrct] = PhantomAnalysis(PhantomEchoesStrct, TransmittedPulsesStruct, PreyFindsStruct, AllParams)
% function [ PhantomAnalysisStrct] = ...
% PhantomAnalysis(PhantomEchoesStrct, TransmittedPulsesStruct, PreyFindsStruct, AllParams)
% 
% the function alayzes the echoes from phnatom reception (conspecifics
% signals echoed from targets and received by the bat
% a phantom will be treateted as a target only if the following conditions
% are followed:
% 1) transmitted signal and received echoe are from the same phase (earch,
% approach, buzz)- by comaring frequencies and durations of the pulses
% (option) - AllP
% 2) the phantom echo is the first echo
% 3) the phantom echo is not masked by 'real echoes' -optional


%%% similiarity thresholds
CheckPhantomSimiliratyFlag = AllParams.BatSonarParams.PhantomEchoesCheckSimiliarity;
freq_th = 4 ; %kHz
duration_th = 1.5e-3/AllParams.SimParams.SampleTime; %1.5msec in samples

%%% the relevant times to analyze (between current own pulse ans next one)
minTime = round(TransmittedPulsesStruct.StartPulseTime + TransmittedPulsesStruct.PulseWidth+1);
maxTime = TransmittedPulsesStruct.StartPulseTime +   ...
    max( TransmittedPulsesStruct.IPItoNextPulse , TransmittedPulsesStruct.PulseWidth+1 );
ReceptionTimes = minTime  : maxTime; 

PhantomTimes = [PhantomEchoesStrct(ReceptionTimes).nTime];
PhantomsIsSimiliar = 0;
IsPhantomDetected = 0;

if ~isempty(PhantomTimes) % is there any relevant Phantom
   
    if CheckPhantomSimiliratyFlag
        %% replace this condition whe it is relevant
        freq_diff = max( ...
            abs( (max([PhantomEchoesStrct(ReceptionTimes).Freqs] - max(TransmittedPulsesStruct.PulseFreqsCommands))) ), ...
            abs( min([PhantomEchoesStrct(ReceptionTimes).Freqs] - min(TransmittedPulsesStruct.PulseFreqsCommands)) ) ); % not good enough
        duration_diff = abs(TransmittedPulsesStruct.PulseWidth - (maxTime-minTime)); % wrong
        if (freq_diff <= freq_th) && (duration_diff <= duration_th) % condition #1
            PhantomsIsSimiliar = 1;
        else % if (freq_diff <= freq_th) & (duration_diff <= duration_th)
            PhantomsIsSimiliar = 0 ;
        end % if (freq_diff <= freq_th) & (duration_diff <= duration_th)
    else % if CheckPhantomSimiliratyFlag
        PhantomsIsSimiliar = 1;
    end % if CheckPhantomSimiliratyFlag
    %%
    if PhantomsIsSimiliar
        PhantomStartTimes = unique([PhantomEchoesStrct(ReceptionTimes).StartnTime]); % sorted by starting times
%         PhantomEndTimes = PhantomEchoesStrct(iST).EndnTime;
        % Condition #2 - the first phantom
        if isempty(PreyFindsStruct.DetectedTimesWithOutInterfernece)
            IsPhantomDetected = 1;
        elseif PhantomStartTimes(1) < min(PreyFindsStruct.DetectedTimesWithOutInterfernece)
            IsPhantomDetected = 1;
        end % if isempty(PreyFindsStruct.DetectedTimesWithOutInterfernece)
    end % if PhantomsIsSimiliar 
end % if ~isempty(PhantomTimes)

if IsPhantomDetected
    %%% dist=(deltaTime*soundspeed/2)
    Est_Dist =  0.5 * AllParams.SimParams.SoundV0 * ...
        (PhantomStartTimes(1)- TransmittedPulsesStruct.StartPulseTime) * AllParams.SimParams.SampleTime;
    % check if teher is more than one phantom at this time and anaylze the
    % storngest
    if numel([PhantomEchoesStrct(PhantomStartTimes(1)).MaxRxPowerDB]) > 1
        [~,curr_max_ind] = max(PhantomEchoesStrct(PhantomStartTimes(1)).MaxRxPowerDB  );
    else
        curr_max_ind=1;
    end % if numel(PhantomEchoesStrct(ReceptionTimes).RxPowersDB > 1)
    estError = 0;
    PhantomAnalysisStrct = struct( ...
        'PhantomIsDetected', 1, ...
        'PhantomPreyNum', PhantomEchoesStrct(PhantomStartTimes(1)).PreyNum(curr_max_ind), ...
        'PhantomRxTime', PhantomEchoesStrct(PhantomStartTimes(1)).nTime(curr_max_ind), ...
        'PhantomEstimatedDistance', Est_Dist, ...
        'PhantomRealDistance', PhantomEchoesStrct(PhantomStartTimes(1)).Distance2PhantomReal(curr_max_ind), ... 
        'PhantomEstimatedAngle', PhantomEchoesStrct(PhantomStartTimes(1)).Angle2Phantom(curr_max_ind) + estError, ...
        'PhantomMaxRxPowerDB',  PhantomEchoesStrct(PhantomStartTimes(1)).MaxRxPowerDB(curr_max_ind), ...
        'PhantomPhaseSignal',PhantomEchoesStrct(PhantomStartTimes(1)).SignalPhase{curr_max_ind} ...
        );
else % if PhantomIsDetected
    PhantomAnalysisStrct = struct([]);
end % if PhantomIsDetected

PhantomAnalysisStrct;