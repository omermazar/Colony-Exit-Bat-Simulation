function [Masking2DetectionStruct] = Interfernce2DetectedPreys( PulsePower,...
                        DetectedPreysVec , CurrEchosFromPreyStruct, InterferenceDirectVec, ...
                        InterferenceFullVec,  InterferenceFullStruct, PulseNum, BAT, BatNumRx, AllParams, FilterBank )
                    
% %  [BAT(BatNum).FindsMaskingStruct ] = Interfernce2DetectedPreys( ...
%                         DetectedPreysVec , CurrEchosFromPreyStruct, InterferenceDirectVec, ...
%                         InterferenceFullVec,  InterferenceFullStruct, AllParams );
%
% Interfernce2DetectedPreys() returns a sturct of the masking to the
% detected prey items including the loclization errors
%  The function is for on-line analysis
% Inputs-   PulsePower - the power of the current pulse [ dB ]
%           DetectedPreysVec - a vector of the detectetd pry items at this
%                   pulse
%           CurrEchosFromPreyStruct - the echos from prey items at this
%               pulse
%           vectors of the direct and full interfernce from other bats (noise level if the is
%               none)
%           InterferenceFullStruct - a struct with whole the times freqs
%               and powers of the interfernces pulses
%           PulseNum - the current tx Pulse the that excited the echoes
%           BAT- struct data of all Bats
%           BatNumRx - the number of the reciveiving bat
%           FilterBank : a struct of the parameters of the filter bank, FilterBank is empty if
%           reciever is not 'FilterBank' 

% Changed: 20Jan2021 - Swarm by omer
%%
%%% START
%%% Parameters setup %%%
% General Params
SampleTime = AllParams.SimParams.SampleTime;
NumberOfPreys = AllParams.SimParams.TotalPreysNumber; 
NumOfSamples = AllParams.SimParams.SimulationTime / SampleTime + 1;
% NumOfTransmittedPulses = BAT.NumOfTimesSonarTransmits;
NoiseLevel = 10.^(AllParams.BatSonarParams.NoiseLeveldB/10);
NoiseLeveldB = AllParams.BatSonarParams.NoiseLeveldB;

% Detection Params
DetectionTH = AllParams.BatSonarParams.PulseDetectionTH;
MinSignal2InterfererTH = AllParams.BatSonarParams.Signal2InterferenceRatio; % dB 
FwdMaskFlag = AllParams.BatSonarParams.MaskingFwdFlag;
FwdMaskingTime = round( AllParams.BatSonarParams.MaskingFwdTime /1000 / SampleTime); % from msec to samples
BckdMaskFlag = AllParams.BatSonarParams.MaskingBckdFlag;
BckdMaskingTime = round( AllParams.BatSonarParams.MaskingBckdTime /1000 / SampleTime); % from msec to samples
MaskingFwdSIRth = AllParams.BatSonarParams.MaskingFwdSIRth; % dB
MaskingBckdSIRth = AllParams.BatSonarParams.MaskingBckdSIRth; % dB
MaskingCorrelationTH = AllParams.BatSonarParams.MaskingCorrelationTH; %dB

MinTimeToDetectEcho = AllParams.BatSonarParams.MinTimeToDetectEcho /1000 / SampleTime; % from msec to samples

% Receiver Type

    ReceiverType = AllParams.BatSonarParams.ReceiverType;
    MaskingCorrelationTH = AllParams.BatSonarParams.MaskingCorrelationTH; %dB
%     CorrelationTimeResolution = AllParams.BatSonarParams.CorrelationTimeResolution; % msec;
    CorrelationTimeResolution = AllParams.BatSonarParams.MaskingFwdTime;
    
    LocalizationErrFlag = AllParams.BatSonarParams.LocalizationErrFlag;  % move to all params ans GUI


%%% Input DATA - 
% interfernce signals
AllInterPowerVec = 10*log10(abs(InterferenceFullVec));
% DirectInterPowerVec = 10*log10(abs(InterferenceDirectVec)); % The vector of the Interferences from other bats' direct calls
% EchosInterPowerVec = 10*log10(abs(InterferenceFullVec - InterferenceDirectVec + NoiseLevel));
% INT = InterferenceFullStruct;

%%% Init Outputs
VecSize = size(DetectedPreysVec);
MaskedPreys = zeros(VecSize);
CurrentPick = zeros(VecSize);
IsPreyMaskedFlag = zeros(VecSize);
IsCoverTimeJamFlag = zeros(VecSize);
CoverTimeJamTimes = zeros(VecSize);
IsFwdMaskingJamFlag = zeros(VecSize);
FwdMaskingJamTimes = zeros(VecSize);
DetecetedPrey2InterferenceRatioDB = zeros(VecSize);
EchoSelfCorrealationMaxdB = -20*ones(VecSize);
EchoMaxInterCorrelationdB = -20*ones(VecSize);
kEchoPrey =0 ;
nSamplesMaxCheck = 100; % [0.5] msec for FsReconstruct=200000
%%% Refernce Outputs Init
RefMaskedPreys = zeros(VecSize);
RefFwdMaskingJamTimes = zeros(VecSize);
RefDetecetedPrey2InterferenceRatioDB = zeros(VecSize);
RefIsFwdMaskingJamFlag = zeros(VecSize);
RefIsPreyMaskedFlag = zeros(VecSize);
RefIsCoverTimeJamFlag = zeros(VecSize);
RefCoverTimeJamTimes = zeros(VecSize);
InterBatNum = 0;
%%
%%% Detected Echos
if  ~strcmp(AllParams.BatSonarParams.ReceiverType, 'FilterBank')
    for nDetPrey = DetectedPreysVec
        kEchoPrey= kEchoPrey+1;
        CurrInd = find(CurrEchosFromPreyStruct.TargetIndex == nDetPrey);
        nEchoTimes = CurrEchosFromPreyStruct.EchoDetailed(CurrInd).EchosnTimesFull;
        nEchoTimes = min(nEchoTimes, NumOfSamples); % protecting the max number of samples
        nEhchoPowers = PulsePower + ...
            10*log10(abs(CurrEchosFromPreyStruct.EchoDetailed(CurrInd).EchoAttenuationFull)) ;
        nEhchoFreqs = CurrEchosFromPreyStruct.EchoDetailed(CurrInd).EchosFreqsFull;
        CurrentPick(kEchoPrey) = max(nEhchoPowers);
        EchoDuration= length(nEchoTimes);
        %%% Check for interfernce for every succesful detection
        
        %%% Reference Check ('MaxPowerTH' - maxPowerDetection)
        %%
        RefMaxTimeCoverInterfernce  = max(AllInterPowerVec(nEchoTimes));
        [RefIsCoverTimeJamFlag(kEchoPrey) , RefCoverTimeJamTimes(kEchoPrey)] = CheckMaskingForMaxPowerTH ...
            (CurrentPick(kEchoPrey), RefMaxTimeCoverInterfernce, MinSignal2InterfererTH, nEchoTimes(1), AllParams);
        % Fwd Masking Interfernce
        FwdMaskStart = max(1,nEchoTimes(1)- FwdMaskingTime) ; % Protect against negative indexes
        [RefMaxFwdMaskingInterfernce, RefFwdSampleNum] = max(AllInterPowerVec( FwdMaskStart : max(nEchoTimes(1)-1, FwdMaskStart) ));
        %     [RefMaxFwdMaskingInterfernce, RefFwdSampleNum] = max(AllInterPowerVec( (nEchoTimes(1)- FwdMaskingTime) : (nEchoTimes(1)-1) ));
        
        [RefIsFwdMaskingJamFlag(kEchoPrey) , RefFwdMaskingJamTimes(kEchoPrey)] = CheckMaskingForMaxPowerTH...
            (CurrentPick(kEchoPrey), RefMaxFwdMaskingInterfernce, MaskingFwdSIRth, nEchoTimes(1), AllParams);
        % check if the Masking is empty and assign Noise level if it is
        if isempty( RefMaxFwdMaskingInterfernce)
            RefMaxFwdMaskingInterfernce = AllParams.BatSonarParams.NoiseLeveldB;
        end % if isempty( RefMaxFwdMaskingInterfernce)
        
        RefIsPreyMaskedFlag(kEchoPrey) = RefIsCoverTimeJamFlag(kEchoPrey) || RefIsFwdMaskingJamFlag(kEchoPrey);
        
        RefDetecetedPrey2InterferenceRatioDB(kEchoPrey) =  CurrentPick(kEchoPrey) - max(RefMaxFwdMaskingInterfernce, RefMaxTimeCoverInterfernce);
        %%
        %%% end of Referece Check
        switch ReceiverType

            %%
            case 'MaxPowerTH'   %  switch ReceiverType
                % This receeveri is the refernce one
                %             CurrentPick(k) = max(nEhchoPowers);
                % Time cover Interference (both echo and direct)
                IsPreyMaskedFlag(kEchoPrey) = RefIsPreyMaskedFlag(kEchoPrey);
                IsCoverTimeJamFlag(kEchoPrey) = RefIsCoverTimeJamFlag(kEchoPrey);
                IsFwdMaskingJamFlag(kEchoPrey)= RefIsFwdMaskingJamFlag(kEchoPrey);
                CoverTimeJamTimes(kEchoPrey) =  RefCoverTimeJamTimes(kEchoPrey);
                FwdMaskingJamTimes(kEchoPrey) = RefFwdMaskingJamTimes(kEchoPrey);
                DetecetedPrey2InterferenceRatioDB(kEchoPrey) = RefDetecetedPrey2InterferenceRatioDB(kEchoPrey);
                InterMinFreq =0;
                InterBW = 0;
                
                %%% end case 'MaxPowerTH'
                %%
            case 'CorrelationDetector'           %  switch ReceiverType
                
                [ SelfCorrelation, SelfCorrTimes, InterefernceCorr, InterCorrTimes,...
                    InterMinFreq, InterBW, InterBatNum] = ...
                    EchosCorrelationsWithTxPulse( PulsePower, CurrEchosFromPreyStruct.EchoDetailed(CurrInd), ...
                    InterferenceFullStruct, FwdMaskingTime,  BAT, BatNumRx, PulseNum, AllParams );
                
                [MaxSelfCorr, MaxSelfIndex] =  max(abs(SelfCorrelation));
                SelfTimeCorr = SelfCorrTimes(MaxSelfIndex);
                [MaxInterCorr, MaxInterIndex] =  max(abs(InterefernceCorr));
                InterTimeCorr = InterCorrTimes(MaxInterIndex);
                InterIndexOfMaxSelfTime = find(InterCorrTimes == SelfTimeCorr);
                if ~isempty( InterIndexOfMaxSelfTime )
                    MaxIndex2Check = min( InterIndexOfMaxSelfTime+nSamplesMaxCheck, length(InterefernceCorr) );
                    MinIndex2Check = max( InterIndexOfMaxSelfTime-nSamplesMaxCheck, 1 );
                    InterfenceCorrWhileSelfMax = max(abs(InterefernceCorr(MinIndex2Check : MaxIndex2Check)));
                else % if ~isempty( InterIndexOfMaxSelfTime )
                    InterfenceCorrWhileSelfMax = NoiseLevel;
                end % if ~isempty( InterIndexOfMaxSelfTime )
                Corr2MaskRatio = 10*log10(max( ...
                    abs(MaxSelfCorr/ MaxInterCorr),...
                    abs(MaxSelfCorr/InterfenceCorrWhileSelfMax) ) ); % PeakDetection Of the self corr
                MinSig2InterferencrRatio = 10*log10(MaxSelfCorr/MaxInterCorr); % the ratio of maxima
                CorrTimeDiff = SelfTimeCorr - InterTimeCorr; % the time diff between the maxima [msec]
                
                [IsCoverTimeJamFlag(kEchoPrey) , CoverTimeJamTimes(kEchoPrey)] = CheckMaskingForCorrelation ...
                    (Corr2MaskRatio, MinSig2InterferencrRatio, CorrTimeDiff, ...
                    MaskingCorrelationTH, CorrelationTimeResolution, nEchoTimes(1));
                
                IsPreyMaskedFlag(kEchoPrey) = IsCoverTimeJamFlag(kEchoPrey);
                DetecetedPrey2InterferenceRatioDB(kEchoPrey) = min(MinSig2InterferencrRatio,Corr2MaskRatio);
                EchoSelfCorrealationMaxdB(kEchoPrey) = 10*log10(abs(MaxSelfCorr));
                EchoMaxInterCorrelationdB(kEchoPrey) = 10*log10(abs(MaxInterCorr));
                
                %             if IsPreyMaskedFlag
                %                 IsCoverTimeJamFlag = 1;
                %                 CoverTimeJamTimes(k) = nEchoTimes(1);
                %             end % if IsPreyMaskedFlag
                FwdMaskingJamTimes(kEchoPrey) = 0; % in correlation there is no meaning for forward masking... We check the time diff
                %             CoverTimeJamTimes(k) = 0;
                
                %%% end case 'CorrelationDetector'
                
                %%
            case 'TimeSequenceTH'  %  switch ReceiverType
                % Init
                FwdMaskingVec = zeros(1,FwdMaskingTime);
                MaskingConv = ones(1,FwdMaskingTime+1);
                
                MinTime2Detect = min(MinTimeToDetectEcho, EchoDuration); % the required num of samples  for detection
                
                % the Detection Window
                
                ks = 1; % index for signal
                
                % init the prameters
                IsDetected = zeros(1, length(nEchoTimes) - MinTime2Detect+1);
                FwdS2I = zeros(1, length(nEchoTimes) - MinTime2Detect+1);
                BwdS2I = zeros(1, length(nEchoTimes) - MinTime2Detect+1);
                for km= nEchoTimes(1:end-MinTime2Detect+1)
                    SignalWindow = [ks: ks+ MinTime2Detect-1];
                    FwdWinnTimes= [max(1,km-FwdMaskingTime):km];
                    CoverWinnTimes = [km: min(km+MinTimeToDetectEcho-1, NumOfSamples) ];
                    % Check For Fwd Masking
                    % Signal2Interference
                    FwdS2I(ks) = max( nEhchoPowers(SignalWindow) ) - ...
                        max(AllInterPowerVec(FwdWinnTimes));
                    BwdS2I(ks) = max(nEhchoPowers(SignalWindow) ) - ...
                        max(AllInterPowerVec(CoverWinnTimes));
                    % Check for fwd Masking
                    IsFwdMaskFlag = FwdS2I(ks) <= MaskingFwdSIRth &...
                        AllInterPowerVec(km) >= NoiseLeveldB + MaskingFwdSIRth;
                    % Check for Bwd Masking
                    IsCoverMaskFlag =  BwdS2I(ks) <= MaskingBckdSIRth & ...
                        AllInterPowerVec(km) >= NoiseLeveldB + MaskingBckdSIRth ;
                    IsDetected(ks) = ~IsFwdMaskFlag & ~IsCoverMaskFlag;
                    
                    ks =ks+1;
                end % for mm
                
                % Results
                DetectedIDx = find(IsDetected);
                IsPreyMaskedFlag(kEchoPrey) = isempty(DetectedIDx);
                IsCoverTimeJamFlag(kEchoPrey) = IsPreyMaskedFlag(kEchoPrey); % No separation between  cover and FwdMasking
                IsFwdMaskingJamFlag(kEchoPrey)= 0; % not relevant
                CoverTimeJamTimes(kEchoPrey) =   IsPreyMaskedFlag(kEchoPrey)* nEchoTimes(1); % zero if no Interfernce
                FwdMaskingJamTimes(kEchoPrey) = 0;
                InterMinFreq =0;
                InterBW = 0;
                %SIR  =
                if IsPreyMaskedFlag(kEchoPrey) % masked
                    DetecetedPrey2InterferenceRatioDB(kEchoPrey) = min([BwdS2I, FwdS2I]);
                else % if IsPreyMaskedFlag(k)
                    DetecetedPrey2InterferenceRatioDB(kEchoPrey) = max([BwdS2I, FwdS2I]);
                end %if IsPreyMaskedFlag(k)
                
                %%% end case 'TimeSequenceTH'
                
                %%
                
                %%% this reciever ignores the masking
                %%
            case 'NoMasking' %
                IsPreyMaskedFlag(kEchoPrey) = 0;
                IsCoverTimeJamFlag(kEchoPrey) = 0; % No separation between  cover and FwdMasking
                IsFwdMaskingJamFlag(kEchoPrey)= 0; % not relevant
                CoverTimeJamTimes(kEchoPrey) =  0; % zero if no Interfernce
                FwdMaskingJamTimes(kEchoPrey) = 0;
                DetecetedPrey2InterferenceRatioDB(kEchoPrey) = CurrentPick(kEchoPrey) - NoiseLeveldB;
                InterMinFreq =0;
                InterBW = 0;
                %%% %%% end case 'NoMasking'
                
                %%
        end  % switch ReceiverType
        
        if IsPreyMaskedFlag(kEchoPrey)
            MaskedPreys(find(~MaskedPreys,1)) = nDetPrey;
        end %if IsPreyMaskedFlag
        
        if RefIsPreyMaskedFlag(kEchoPrey)
            RefMaskedPreys(find(~RefMaskedPreys,1)) = nDetPrey;
        end %if IsPreyMaskedFlag
        
        
    end % for nDetPrey = DetectedPreysVec
else % if  ~strcmp(AllParams.BatSonarParams.ReceiverType, 'FilterBank')

end % if  ~strcmp(AllParams.BatSonarParams.ReceiverType, 'FilterBank')
%%
%%% Calculate lcalization error (DF, Distance, RElative pos)
% CalculateLocalizationErros function nedees; the Rxpower, SIR, Angle,
% Distance, TransmittedPulseTime

if LocalizationErrFlag
    [DFErr, RangeErr, RelativeDirectionErr] = CalculateLocalizationErrors(...
        BAT(BatNumRx), CurrEchosFromPreyStruct, DetectedPreysVec, ...
        DetecetedPrey2InterferenceRatioDB, CurrentPick, PulseNum, AllParams);
else % if LocalizationErrFlag
    DFErr =0;
    RangeErr = 0;
    RelativeDirectionErr = 0;
end % if LocalizationErrFlag
%%
%%% OUTPUT
TotalMaskingTimes = nonzeros(unique([CoverTimeJamTimes, FwdMaskingJamTimes]) )';
RefTotalMaskingTimes = nonzeros(unique([RefCoverTimeJamTimes, RefFwdMaskingJamTimes]) )';

IsAnyPreyMasked = (sum(IsPreyMaskedFlag) > 0 );

% add Vectors of Acoustic Classifier - New 15Jun22 - Omer
classifier_missed_targets    = [];
classifier_jammed_targets    = [];
direct_jammed_targets        = [];
classifier_jammed_times      = [];
FalseAlarmsClassified_nTimes = [];
FalseAlarmsDistances         = [];
FalseAlarmsDoa               = [];
classifierResultsAll         = struct();
DetectedCalissified_nTimes   = [];
classifier_jammed_pks        = [];

% The Output
% NEW 08Dec2022 - Choose only unmasked targets
ixMasked = ismember(DetectedPreysVec, nonzeros(MaskedPreys)'); 
RxDetectedOnly = CurrentPick(~ixMasked);

Masking2DetectionStruct = struct(...
        'PulseNum', PulseNum, ...
        'DetectedPreys', DetectedPreysVec,...
        'DetectedPreysRxPower', RxDetectedOnly, ...% CurrentPick, ...
        'IsAnyPreyMasked', IsAnyPreyMasked, ...
        'MaskedPreys', nonzeros(MaskedPreys)', ...
        'masked_prey_idx',[] , ...
        'TotalMaskingTimes', TotalMaskingTimes,...
        'IsCoverTimeJamFlag' , IsCoverTimeJamFlag,...
        'CoverTimeJamTimes' , CoverTimeJamTimes, ...
        'IsFwdMaskingJamFlag', IsFwdMaskingJamFlag,...
        'FwdMaskingJamTimes', FwdMaskingJamTimes ,...
        'DetecetedPrey2InterferenceRatioDB', DetecetedPrey2InterferenceRatioDB, ...
        'SelfCorrealationMaxdB', EchoSelfCorrealationMaxdB, ...
        'InterCorrelationMaxdB', EchoMaxInterCorrelationdB, ...
        'Ref_MaskedPreys', nonzeros(RefMaskedPreys)',  ...
        'Ref_TotalMaskingTimes', RefTotalMaskingTimes, ...
        'Ref_FwdMaskingJamTimes', RefFwdMaskingJamTimes, ...
        'Ref_DetecetedPrey2InterferenceRatioDB', RefDetecetedPrey2InterferenceRatioDB,...
        'DetectionsDFerror', DFErr,...
        'DetectionsRangeErr', RangeErr,...
        'DetectionsRelativeDirectionErr', RelativeDirectionErr, ...
        'InterMinFreq', InterMinFreq, ...
        'InterBW',   InterBW, ...
        'InterBatNum', InterBatNum, ...
        'FB_unmasked_prey', [], ...
        'FB_unmasked_delays', [], ...
        'FB_unmasked_times',[], ...
        'FB_detected_masking_delays', [], ...
        'FB_detected_masking_powers', [] , ...
        'FB_detected_Prey_times', [] , ...
        'FB_estimated_masker', [], ...
        'ClassifierMissedTargets', classifier_missed_targets, ...
        'ClassifierMaskedTargets', classifier_jammed_targets , ...
        'DetectionMaskedTargets', direct_jammed_targets, ...
        'ClassifierMaskedTimes', classifier_jammed_times ,...
        'ClassifierFalseAlarmsTimes', FalseAlarmsClassified_nTimes ,...
        'ClassifierFalseAlarmsDistances', FalseAlarmsDistances ,...
        'ClassifierFalseAlarmsDOA', FalseAlarmsDoa ,...
        'ClassifierDetecdtedAndTRUE', DetectedCalissified_nTimes, ...
        'ClassifierFullResults',  classifierResultsAll ...   
        );
    
    
    
end % function

%%% End of main function %%%%%

%%
function [IsMaskedFlag, MaskingTime] = ...
    CheckMaskingForMaxPowerTH (MaxSignalPower, MaxInterferencePower, SIRTH, EchoStartTime, AllParams)

NoiseLeveldB = AllParams.BatSonarParams.NoiseLeveldB;

IsMaskedFlag = 0;
MaskingTime = 0;
% if the interferfernce is weaker than the noise, the TH is set do
%     detection lenvel (SNR)
if MaxInterferencePower - NoiseLeveldB <= 1
    SIRTH = min(SIRTH, AllParams.BatSonarParams.PulseDetectionTH); % if it's fwd time
    % else % if MaxInterferencePower - NoiseLeveldB <= 1
    
end % if MaxInterferencePower - NoiseLeveldB <= 1

if MaxSignalPower <= MaxInterferencePower + SIRTH
    IsMaskedFlag = 1;
    MaskingTime = EchoStartTime;
end % if MaxSignalPower

end % function  CheckMaskingForMaxPowerTH
%%

%%
function [IsCorrJamFlag , CorrJamTimes] = CheckMaskingForCorrelation ...
    (Corr2MaskRatio, MinSig2InterferencrRatio, CorrTimeDiff, ...
    MaskingCorrelationTH, CorrelationTimeResolutionTH, EchoStartTime)
% reurnns a flag whether a jam happens, and the time of the jamming (the
% start echo pulse

IsCorrJamFlag = 0;
CorrJamTimes = 0;
 % check masking from the 'skirt of the correlation'
if Corr2MaskRatio <= MaskingCorrelationTH
    IsCorrJamFlag = 1;
    CorrJamTimes = EchoStartTime;
    
% check  masking between ''
elseif (MinSig2InterferencrRatio <= MaskingCorrelationTH) && ...
        (CorrTimeDiff <= CorrelationTimeResolutionTH)
    IsCorrJamFlag = 1;
    CorrJamTimes = EchoStartTime;
end % if Corr2MaskRatio <= MaskingCorrelationTH

if IsCorrJamFlag 
    CorrJamTimes = EchoStartTime;   
end % if IsCorrJamFlag 
end % function CheckMaskingForCorrelation

%%
