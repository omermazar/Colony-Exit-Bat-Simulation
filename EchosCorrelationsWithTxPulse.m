function [ SelfCorrelation, SelfCorrTimes, InterferenceCorrelation, InterferenceCorrTimes, ...
    InterMinFreq, InterBW, InterBatNum ] = ...
    EchosCorrelationsWithTxPulse( PulsePowerdB, EchoDetailedStruct, ...
                InterferenceFullStruct, FwdMaskingTime,  BAT, BatNumRx, OwnPulseNum, AllParams )
%EchosCorrelationsWithTxPulse
%   the function reconstruct the pulse, the ecgos and the interdrence
%   then calculttes the cros correlation between the echo and the pulse
%   and betweent the Interfrence and thr pulse

%%% Re- Contructing the echo and the Tx pulse
%%
SampleTime = AllParams.SimParams.SampleTime;
NumOfSamples = AllParams.SimParams.SimulationTime / SampleTime + 1;

FsReconstruct = 200*1e3; % the sampling rate of the reconstruced signals
DeltaFreq = 5; % To Be Changed

PulsePower= 10^(PulsePowerdB/10);
NoiseLevel = 10^(AllParams.BatSonarParams.NoiseLeveldB/10);
OwnPulseMinMaxFreqs = BAT(BatNumRx).TransmittedPulsesStruct(OwnPulseNum).ChirpMinMaxFreqs;
EchoMinFreq = OwnPulseMinMaxFreqs(1);
EchoMaxFreq = OwnPulseMinMaxFreqs(2);

EchoTime = [EchoDetailedStruct.EchosnTimesFull];
% EchoDuration = length(EchoDetailedStruct.EchosnTimesFull);
% EchoMinFreq =  min(EchoDetailedStruct.EchosFreqsFull);
EchoAtten =  max(abs(EchoDetailedStruct.EchoAttenuationFull)); % later - use interpulation
% if EchoDuration>1
%     EchoWidth = EchoDuration * SampleTime; %msec
%     EchoMaxFreq = max(EchoDetailedStruct.EchosFreqsFull);
% else % if PulseDuration>1
%    EchoWidth = AllParams.BatSonarParams.PulseTimeShort ; %msec
%    EchoMaxFreq =  EchoMinFreq + DeltaFreq ; 
% end %if PulseDuration>1
EchoWidth = BAT(BatNumRx).TransmittedPulsesStruct(OwnPulseNum).PulseWidth*SampleTime;

t= 0 : 1/FsReconstruct : EchoWidth;
%%% Reconstructin the tx pulse
% The parameters are the same as the echo
TxPulseChirp =  chirp(t, EchoMaxFreq*1000, max(t), EchoMinFreq*1000,'logarithmic');
EchoChirp = TxPulseChirp * PulsePower * EchoAtten;

InterMinFreq = 0;
InterBW = 0;
InterBatNum = 0; 

%%% TEST
% EchoChirp = chirp(t,70000, max(t), 40000,'logarithmic');
%%% spectrogram
%spectrogram(EchoChirp,64,63,64,FsReconstruct,'yaxis')

%%
%%% The Self Correlation  - done at the end
%%
% [SelfCorrelation,SelfLags] = xcorr(EchoChirp, TxPulseChirp);
% IntegratioTime = 0.1*1e-3; % 100 micro
% Kcorr = IntegratioTime * FsReconstruct;
% SelfCorrelation = SelfCorrelation/Kcorr;
% SelfCorrTimes = min(EchoTime) + SelfLags/ FsReconstruct/ SampleTime;
%%

%%% Reconstruct the Masking Signal
%%
RelInterferenceTime = max(min(EchoTime)-FwdMaskingTime, 1) : min(max(EchoTime),NumOfSamples);
CurrInterStruct = InterferenceFullStruct(RelInterferenceTime);
InterferenceTimes = nonzeros([CurrInterStruct.Times]);
% Check if there is interfernce signal during Times, if there is-
% reconstruct
if ~isempty(InterferenceTimes)
    MaxSampleInterPower = zeros(size(RelInterferenceTime));
    MaxnInd = zeros(size(RelInterferenceTime));
    
    %%% Reconstructing the chirp with the maximum Power
    for nn = 1:length(RelInterferenceTime)
        if ~isempty([CurrInterStruct(nn).Power])
            [MaxSampleInterPower(nn), MaxnInd(nn)] = max( abs([CurrInterStruct(nn).Power]) );
        end % if ~isempty(InterferenceTimes(nn))
    end % for nn
    [MaxInterPower, maxPowerInd] = max(MaxSampleInterPower);
    InterBatNum = CurrInterStruct(maxPowerInd).BATTxNum( MaxnInd(maxPowerInd) );
    InterPulseNum = CurrInterStruct(maxPowerInd).BatTxPulseNum( MaxnInd(maxPowerInd) );
    
    InterferenceFreqs = BAT(InterBatNum).TransmittedPulsesStruct(InterPulseNum).ChirpMinMaxFreqs; 
    InterMaxFreq = max(InterferenceFreqs);
    InterMinFreq = min(InterferenceFreqs);
   % BatVariance = 0.05*InterMinFreq; % difference between bats
%     InterMaxFreq = InterMaxFreq; %+ BatVariance*randn(1);
%     InterMinFreq = InterMinFreq; %+ BatVariance*randn(1);
    InterBW = BAT(InterBatNum).TransmittedPulsesStruct(InterPulseNum).PulseWidth;
%     InterBW = InterMaxFreq- InterMinFreq;
    InterStartTime = RelInterferenceTime(maxPowerInd);
    InterTimeWidth = InterBW*SampleTime; %[msec]
%     if InterBW > 0 
%         InterTimeWidth = abs(minInd - maxInd) * SampleTime;
%     else % if InterBW > 0
%         InterTimeWidth = AllParams.BatSonarParams.PulseTimeShort;
%         InterMaxFreq =  InterMinFreq + DeltaFreq ; 
%     end % if InterBW > 0
    ti = 0 : 1/FsReconstruct : InterTimeWidth;
%     SampleInterPower = mean(MaxSampleInterPower);
    InterChirp = MaxInterPower * chirp(ti, InterMaxFreq*1000, max(ti), InterMinFreq*1000,'logarithmic');

% if there is no-interfernce, the chirp is noise
else % if ~isempty(InterferenceTimes
    InterChirp = NoiseLevel*randn(size(TxPulseChirp));
    InterStartTime = min(EchoTime); 
end % if ~isempty(InterferenceTimes)


%%% The correlations
%%
%%% pading the short chirp
LTx= length(TxPulseChirp);
LInt = length(InterChirp);
PadSize = abs(round((LTx- LInt)/2));
if LTx > LInt
    InterChirp = padarray(InterChirp', PadSize, 0,'both')';
elseif LTx < LInt % if LTx > LInt
    EchoChirp = padarray(EchoChirp', PadSize, 0,'both')';
    TxPulseChirp = padarray(TxPulseChirp', PadSize, 0,'both')';
end % if LTx > LInt

%%% The Self Correlation
[SelfCorrelation,SelfLags] = xcorr(EchoChirp, TxPulseChirp);
% try % add to AllParams
    % AllParams.BatSonarParams.CorrelationTimeResolution is the integration time in msec
    IntegrationTime = AllParams.BatSonarParams.CorrelationTimeResolution*1e-3; %  / AllParams.SimParams.SampleTime;
% catch
%    IntegrationTime = 0.1*1e-3; % 100 micro
% end% try
% IntegrationTime = 0.1*1e-3; % 100 micro
Kcorr = IntegrationTime * FsReconstruct;
SelfCorrelation = SelfCorrelation/Kcorr;
SelfCorrTimes = min(EchoTime) + SelfLags/ FsReconstruct/ SampleTime;

%%% Clculate Cross Correlation
[InterferenceCorrelation, CrossLags] = xcorr(InterChirp, TxPulseChirp);
InterferenceCorrelation = InterferenceCorrelation/Kcorr;
InterferenceCorrTimes = InterStartTime + CrossLags/ FsReconstruct/ SampleTime;
%%

%%% Sum Vector
%%
% minTime= min( min(CrossCorrTimes), min(SelfCorrTimes) );
% maxTime= max( max(CrossCorrTimes), max(SelfCorrTimes) );
% SumTimes = minTime: SampleTime: maxTime;
% SumCorr = NoiseLevel*ones(size(SumTimes));

%%
end

