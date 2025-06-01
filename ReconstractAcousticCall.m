function [chirpCall] = ReconstractAcousticCall(TransmittedPulsesStruct, fs, nSamples)

% chirpCall = ReconstractAcousticCall(TransmittedPulsesStruct)
%
% returns the call of the bat as acousic signal sampled at fs, 
% TransmittedPulsesStruct is the struct:
% BAT(BatNum).TransmittedPulsesStruct(PulseNum) with the fields: 
%               PulseNum: 2
%             PulsePower: 110
%         StartPulseTime: 323
%             PulseWidth: 14
%       ChirpMinMaxFreqs: [34.3316 42.3316]
%     PulseFreqsCommands: [41.7030 41.0836 40.4735 39.8725 39.2803 38.6970 38.1223 37.5562 36.9985 36.4490 35.9077 35.3745 34.8491 34.3316]
%         IPItoNextPulse: 200
%          PulseDuration: 14
%      ManueverTypePower: 'HuntingRegularHunt'
%          JARPulsesLeft: 0
% SimSampleTime: Params.SimParams.SampleTime

%% Input
% if nargin <3
%     Sim_SampleTime = 0.5e-3; 
% end
% if nargin <2
%     fs = 200e3; 
% end

%% Main
MaxFreq = TransmittedPulsesStruct.ChirpMinMaxFreqs(2);
MinFreq = TransmittedPulsesStruct.ChirpMinMaxFreqs(1);
%EchoDur = TransmittedPulsesStruct.PulseWidth * Sim_SampleTime;
PulsePower = 10^(TransmittedPulsesStruct.PulsePower/10);

%t = 0 : 1/fs : EchoDur;
t = 0 : 1/fs : 1/fs*(nSamples-1);
% The chirp
chirpCall = sqrt(PulsePower) * chirp(t, MaxFreq*1000, max(t), MinFreq*1000, 'logarithmic');

%%% add Fade-in fade=out
fadeFlag = 1;

if fadeFlag

    % Parmeters
    fade_percentage = 10; % Percentage of pulse width for fade-in and fade-out
    fade_samples = floor(fade_percentage / 100 * nSamples); % Number of samples for fade-in and fade-out
    fadeMethod = 'Linear'; % 'Gaussian'; % 'Linear'
    
    switch fadeMethod
        case 'Gaussian'
            % Create Gaussian windows for fade-in and fade-out
            gaussian_window = gausswin(2*fade_samples + 1);
            fade_in  = gaussian_window(1:fade_samples);
            fade_out = gaussian_window(fade_samples+2:end);
        case 'Linear'
            fade_in  = linspace(0, 1, fade_samples)';
            fade_out = linspace(1, 0, fade_samples)';
    end % switch

    chirpCall(1:fade_samples) = chirpCall(1:fade_samples) .* fade_in';
    chirpCall(end-fade_samples+1:end) = chirpCall(end-fade_samples+1:end) .* fade_out';
end % if fadeFlag


