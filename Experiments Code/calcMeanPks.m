function [pks_mean, pks_std, pks, pks_times, echoes_mean, echoes_std, echo_pks, echo_times] = calcMeanPks(sig, fs, echo_flag)

%  [pks_mean, pks_std, pks, pks_times] = calcMeanPks(sig, fs, echoFlag)
%   Returnns the mand and std of detectected calls is sig. also returns the
%   peaks' amplitude and time.
% if echoFlag tre - returns the means of the echoes after the calls
% input: 
% path_name = 'C:\audioRecordings\RM_target_strength\results';
% file_name = 'ref_100cm_inFront009.wav'; % for the refencence calculation
% file_name = 'ts_100cm_rotate013.wav'; % for the rotating target echo 
% file_name = 'ts_100cm_static011.wav'; % for the stattic target echo 
%
% [sig, fs] = audioread(fullfile(path_name, file_name));
%% BPF
f_low  = 25e3; % 
f_high = 45e3; % f_high = 49500;

nOrder = 8;
bpf_filter = designfilt('bandpassiir', 'FilterOrder',nOrder, ...
    'HalfPowerFrequency1', f_low, 'HalfPowerFrequency2', f_high, 'SampleRate', fs); 
sig_filt = filtfilt(bpf_filter, sig); 

% time vec
t = [1:numel(sig)]/fs;

%% calc the envelope
env_sig = envelope(sig_filt);
env_sig = smooth(env_sig, 15);

%% peak detection 
% realtive to noise
% prc_peak = 90; 
% detect_TH = prctile(env_sig, prc_peak) * sqrt(1000);
% manually set
% detect_TH = 0.25; % for mic infornt the speaker
detect_TH = 0.05; % for mic and speakear inline (for echoes) echoes

% minimum IPI between detection
if echo_flag
    min_IPI_TH = 4e-3*fs;
else
    min_IPI_TH = 50e-3*fs;
end % if 

% the prominance
min_pr     = 0.5 * detect_TH;

%%%% find peaks with the relvant setection TH
[pks, loc, widths, pr] = findpeaks(env_sig, 'MinPeakHeight',  detect_TH, 'MinPeakProminence', min_pr, ...
'MinPeakDistance',min_IPI_TH);
% the detection times
pks_times = t(loc);

% calc the mean and std
pks_mean = mean(pks);
pks_std  = std(pks);


figure
hold on; grid on;
plot(t, sig, 'b', DisplayName='input')
plot(t, env_sig, 'k', DisplayName='envelope')
plot(pks_times, pks, 'og', DisplayName='pks')

%% find the echo of each call 
if ~echo_flag
    return
end % if echo flag

% params
expected_echo_times = [2,6]*1e-3;
echo_th = 0.005; 
echo_dist = 4e-3*fs; % 4 ms
n_calls = numel(pks);

% init
echo_pks   = nan(1,n_calls);
echo_times = nan(1,n_calls);

for k = 1:n_calls
    ix_echo = t > pks_times(k) + expected_echo_times(1) & ...
              t < pks_times(k) + expected_echo_times(2);
    curr_env = env_sig(ix_echo);

    [p, loc] =max(curr_env);
    echo_pks(k) = p;
    echo_times(k) = pks_times(k) + expected_echo_times(1) + loc/fs;    
end

% calc the mean and std
echoes_mean = mean(echo_pks);
echoes_std  = std(echo_pks);

plot(echo_times, echo_pks, '*r', DisplayName='echoes')




