function [NumOfPreyDetected, DetectedPreysVec, DetectedPreyTimes, ClutteredPreyVec, ClutteredTimes] = ...
    IsPreyDetected(ObsPreyFlag, EchosPreyStruct, TransmittedPulseStruct, ...
    NumObsDetected, CurrEchoesFromObsStruct, AllParams)

% function [NumOfPreyDetected, DetectedPreysVec, DetectedPreyTimes] = ...
%     IsPreyDetected(EchosPreyStruct, TransPulsePower, CurrentIPI, AllParams)

% this function will check whether and which prey is detected
% Detection Conditions:
%   1) EchosPower > Detection Thresh Hold
%   2) Delay Time <= Current IPI (without Ambiguty of types 1,2)
% Outputs:
%   NumOfPreyDetected - the total number of the prey that wewre detected
%   DetectedPreysVec - vectorof the preys that were detectd
% Inputs:
%   EchosPreyStruct - the struct of the echos from all thr preys
%   PulseDetectionTH - the minimpum detcted pulse enrgy (dBm)
%   IsObsDetected - A flag indicatin that an obstacle is detected in the
%           current pulse
%   CurrEchoesFromObsStruct - the last obstacles detection details
%
%%%  Change: 26/02/20 - adding clutter mode %%%
%%% Change: Jan 2021 - addinng ConspecificDetection
%%% Cahnge: May 2022 - Fixing ClutterMode for Acoustics

PulseDetectionTH = AllParams.BatSonarParams.PulseDetectionTH;
NoiseLeveldB = AllParams.BatSonarParams.NoiseLeveldB;
SampleTime = AllParams.SimParams.SampleTime;

NumOfPreyDetected = 0;
DetectedPreysVec = zeros(1,EchosPreyStruct.NumOfEchos);
DetectedPreyTimes = zeros(1,EchosPreyStruct.NumOfEchos);
NumOfPreyCluttered = 0;
ClutteredPreyVec = zeros(1,EchosPreyStruct.NumOfEchos);
ClutteredTimes = zeros(1,EchosPreyStruct.NumOfEchos);

TransPulsePower = TransmittedPulseStruct.PulsePower;
CurrentIPI =  TransmittedPulseStruct.IPItoNextPulse;
TxDur =  TransmittedPulseStruct.PulseDuration * SampleTime;
TerminalFreq =  TransmittedPulseStruct.ChirpMinMaxFreqs(1);
TxBW = TransmittedPulseStruct.ChirpMinMaxFreqs(2) - TerminalFreq;
CorrelationTimeResolution = AllParams.BatSonarParams.CorrelationTimeResolution;

%% New: detect by corr  %%% Change May2022
DetectByCorrFlag = 0;
CorrGainRef      = 0;
CorrGaindB       = 0;
if ~strcmp(ObsPreyFlag, 'Obs')
    %%% Prey Only
    if sum(strcmp('DetectionByCorr', fields(AllParams.BatSonarParams)))
        if  AllParams.BatSonarParams.DetectionByCorr
            DetectByCorrFlag = 1;
            CorrGainRef = AllParams.BatSonarParams.CorrGainRefdB; % 15dB auto-corr of search calls with duration 7ms
            % calculate autocorrelation
            [CorrGaindB] = MyAutoCorr(TerminalFreq, TxBW, TxDur,  CorrelationTimeResolution);

        end % if  AllParams.BatSonarParams.DetectionByCorr
    end % if sum(strcmp())
end % if strcmp(ObsPreyFlag, 'Prey')

%% Clutter Interference
IsObsDetected = NumObsDetected > 0;

if ~strcmp(ObsPreyFlag, 'Obs') && ~AllParams.SimParams.AcousticsCalcFlag
    if sum(strcmp('ClutterMode', fields(AllParams.BatSonarParams)))

        if  AllParams.BatSonarParams.ClutterMode && IsObsDetected
            % verify that the clutter is from the same pulse
            same_pulse_flag = ...
                CurrEchoesFromObsStruct.TransmittedPulseTime == EchosPreyStruct.TransmittedPulseTime;
            if same_pulse_flag
                c_atten_all = [CurrEchoesFromObsStruct.EchoDetailed.EchoAttenuationFull];
                c_times_all = [CurrEchoesFromObsStruct.EchoDetailed.EchosnTimesFull];
                % sum signals by rx times
                c_times = unique(c_times_all);
                c_atten = nonzeros(accumarray(c_times_all'-min(c_times_all)+1, c_atten_all,[], @rssq));
                clutter_powers = TransPulsePower + 10*log10(c_atten);
            else % if same_pulse_flag
                clutter_powers = NoiseLeveldB-10;
            end % if same_pulse_flag
        end % if AllParams.BatSonarParams.ClutterMode
    end %  if sum(strcmp('DetectionByCorr', fields(AllParams.BatSonarParams)))
end % if strcmp(ObsPreyFlag, 'Prey')

%%
% the Echos in EchosPreyStruct are orderd by time not by prey
for nPrey= 1:EchosPreyStruct.NumOfEchos
    isDetected = 0;
    [maxPower, maxIdx] = max(abs(EchosPreyStruct.EchoDetailed(nPrey).EchoAttenuationFull)); % the peak

    %     nPulsePower = TransPulsePower + 10*log10(abs(max(EchosPreyStruct.EchoDetailed(nPrey).EchoAttenuationFull)));
    nPulsePower = TransPulsePower + 10*log10(maxPower);
    nPulseTime = EchosPreyStruct.EchoDetailed(nPrey).EchosnTimesFull(maxIdx);

    %%%% Change May2022
    switch ObsPreyFlag
        case {'Prey', 'Consp'}
            %% Prey Detecetions only in the current IPI ( Condtion #2)
            if EchosPreyStruct.EchosTimes(nPrey) <= CurrentIPI % Condtion #2

                if nPulsePower+CorrGaindB >= NoiseLeveldB + PulseDetectionTH + CorrGainRef % Condition #1
                    isDetected = 1;
                end % if

                %% clutter check
                %%%% Change May2022
                if isDetected && AllParams.BatSonarParams.ClutterMode && ~AllParams.SimParams.AcousticsCalcFlag

                    % the relevant making window around the peak:
                    MaskingStartTime = nPulseTime - AllParams.BatSonarParams.MaskingFwdTime*1e-3/SampleTime;
                    MaskingEndTime = nPulseTime + AllParams.BatSonarParams.MaskingBckdTime*1e-3/SampleTime;
                    idx_clutter_masking = [c_times > MaskingStartTime] & [c_times < MaskingEndTime];

                    % the maximum clutter iin the window:
                    if any(idx_clutter_masking)
                        max_clutter = max(clutter_powers(idx_clutter_masking));
                    else %if ~isemapty(idx_clutter_masking)
                        max_clutter = NoiseLeveldB-10;
                    end % if ~isemapty(idx_clutter_masking)

                    % checking peaks:
                    if nPulsePower < max_clutter+PulseDetectionTH % XXX %
                        isDetected = 0;
                        clutter_flag = 1;
                        NumOfPreyCluttered = NumOfPreyCluttered+1;
                        ClutteredPreyVec(NumOfPreyCluttered) = EchosPreyStruct.TargetIndex(nPrey) ; % The prey ID that was detecet
                        ClutteredTimes(NumOfPreyCluttered) = nPulseTime;
                    end % nPulsePower < max_clutter+PulseDetectionTH


                end %if isDetected && AllParams.BatSonarParams.ClutterMode
            end %if EchosPreyStruct.EchosTimes(nPrey) <= CurrentIPI

        case 'Obs'

            if nPulsePower+CorrGaindB >= NoiseLeveldB + PulseDetectionTH + CorrGainRef % Condition #1
                isDetected = 1;
            end % if

    end % switch ObsPreyFlag

    %%
    if isDetected
        NumOfPreyDetected = NumOfPreyDetected +1;
        DetectedPreysVec(NumOfPreyDetected) =  EchosPreyStruct.TargetIndex(nPrey) ; % The prey ID that was deteceted
        DetectedPreyTimes(NumOfPreyDetected) = round(EchosPreyStruct.EchosTimes(nPrey)) + EchosPreyStruct.TransmittedPulseTime;
    end % if isDetected

end % for nPrey= 1:EchosPreyStruct.NumOfEchos

%% OutPut
DetectedPreysVec  = nonzeros(DetectedPreysVec)';
DetectedPreyTimes = nonzeros(DetectedPreyTimes)';
ClutteredPreyVec  = nonzeros(ClutteredPreyVec)';
ClutteredTimes    = nonzeros(ClutteredTimes)';
