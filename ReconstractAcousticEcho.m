function [AcousticEcho] = ReconstractAcousticEcho(EchoDetailed, chirpCall, fs, nSamples)

% [AcousticEcho] = ReconstractAcousticEcho(EchoDetailed, chirpCall)
%
% returns the echo as acousic signal sampled at fs, 
% EchoDetailed is the struct:
% BAT(BatNum).EchosFromObsStruct(PulseNum).EchoDetailed(k) with the fields: 
%       TargetIndex: 1
%       EchoAttenuationFull: [1×PulseWidth double]
%       EchosnTimesFull: [358 359 360 361 362 363 364 365 366 367 368 369 370 371]
%       EchosFreqsFull: [41.7030 41.0836 40.4735 39.8725 39.2803 38.6970 38.1223 37.5562 36.9985 36.4490 35.9077 35.3745 34.8491 34.3316]
% 
% SimSampleTime: Params.SimParams.SampleTime

%% NOV2021

dur = numel(EchoDetailed.EchosnTimesFull);
t0 = linspace(EchoDetailed.EchosnTimesFull(1), EchoDetailed.EchosnTimesFull(end)+1, dur);
if dur == 1
    t0(2) = t0+fs*nSamples;
    EchoDetailed.EchoAttenuationFull = EchoDetailed.EchoAttenuationFull*ones(1,2);
end
%  tq = EchoDetailed.EchosnTimesFull(1) + (0 : 1/fs : 1/fs*(nSamples-1) );
if ~isinf(t0(1))
    tq = linspace(t0(1), t0(end), nSamples);
    AttenMask = interp1(t0, EchoDetailed.EchoAttenuationFull, tq); % Power (signal level)
    AttenMask = sqrt(AttenMask); % signal
else
    AttenMask = zeros(size(chirpCall));
end % if ~isinf(t0(1))

AcousticEcho = chirpCall.*AttenMask;

% figure
% hold on
% plot(t0, EchoDetailed.EchoAttenuationFull, 'or')
% plot(tq, AttenMask, '.-')
% plot(tq, AcousticEcho)