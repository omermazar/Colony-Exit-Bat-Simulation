function [CorrGaindB] = MyAutoCorr(TerminalFreq, TxBW, TxDur,  CorrelationTimeResolution)
% inputs:   TerminalFreq- thr minimum freq of the pulse [kHz]
%           TxBW - the band-width [kHz]
%           TxDur - the duration of the pulse [s]
%           Fs - sampling (Reconstruct) frequency for correlation calculation (Fs = 200e3;
%           CorrelationTimeResolution - 0.1ms = AllParams.BatSonarParams.CorrelationTimeResolution
Fs = 20e3;
MaxFreq = TerminalFreq + TxBW; %[kHz]
tChirp=  0 : 1/Fs : TxDur; % sec

TxChirp = 1* chirp(tChirp, MaxFreq*1000, ...
    max(tChirp), TerminalFreq*1000,'logarithmic');

AutoCorrelation = xcorr(TxChirp, TxChirp);
IntegrationTime = CorrelationTimeResolution*1e-3; %  / AllParams.SimParams.SampleTime;
Kcorr = IntegrationTime * Fs;
CorrGaindB = 10*log10(max(abs(AutoCorrelation))/Kcorr);
