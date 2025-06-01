function [AllFiltersParams]  = ...
    FilterBank_Acous_Resp(num_of_filters, fs, f_low, f_high, response_time, fft_resolution)

% this function calucates the impulse respnse of gammtone filter bank:
%     gt = a*t.^(n-1) .* exp(-2*pi*b*t) .* cos(2*pi*fc*t);

% output: AllFiltersParams a struct with the fields
%   filter_fc - the center freqiencies of the filters
%   filter_timevec - vector of times of the impulse respnse ([0:ts:imp_toltal_time])
%   filter_impulse_resp  = matrix of the impulse respnse of each filter in the bank
%                   impulse_resp(k,i) is the impluse response of fc(i), in time t(i)
%   filter_cross_talk = matrix of the fft's of each imulse response in all
%   filter_delay_times - vecor of the the delay times of the maxumum impulse resopnse for each fc (samples)
%                   frequencies of the filter bank: filter_cross_talk(k,i)
%                   is the the fft value of fc(k) in the the freq fc(i)
%   fft_freqs = vector of the frequencies input for the fft: [f_low:fft_resolution:f_high] Hz
% inputs:
%   num_of_filters = number of filters in the bank 
%   fs = Sampling frequency (default: 250e3, ts=4 micro seconds)
%   f_low = lowest freq in the filterbank (center freq) (default: 10e3)(Hz)
%   f_high - highest freq in the filter bank (default(80e3) (Hz)
%   response_time - the maximum time of the impluse response (500e-6)
%   fft_resolution - the resolution of the fft output (1e3)
% Default Params:

%%% Nov2021
% 80 chaneles (5 - 100 kHz): fc = 5702*2.^(k_filt/k_norm); freqs 6 - 100 kHz;
% 1) BPF: B = 0.15*fc, Gammatone impulse response: gt = a*t.^(n_order-1) .* exp(-2*pi*B*t) .* cos(2*pi*fc(k)*t); Order -4 
% 2) Rectifier: y=x. if x>0, y=0 if x<0 *
% 3) LPF: Freq = 8khz; Order = 4, Butterworth *
% 4) Tao = delay shift by the maximum peak of each filter *
% 5) Summation *
% * steps 2-5 are applied in the detector as the Retifier is not linear

% REF:
%   Sanderson, M. I., Neretti, N., Intrator, N. & Simmons, J. A. Evaluation of an auditory 936 model for echo delay accuracy in wideband biosonar. J. Acoust. Soc. Am. 114, 1648–937 1659 (2003)
%   Boonman, A. & Ostwald, J. A modeling approach to explain pulse design in bats. Biol. 934 Cybern. 97, 159–172 (2007).

%% general parameters
% num_of_filters = 80;
% fs = 200e3;
% f_low = 10e3;
% f_high = 40e3;
% response_time = 0.5e-3; % good for f_low = 20kHz
% fft_resolution = 1e3;

ts = 1/fs; 
k_filt = 1:num_of_filters;
k_norm = num_of_filters/ log2(f_high/5703);
fc = 5702*2.^(k_filt/k_norm); % THe cEntral Frewquncy of each BPF, total of num_of_filters between ~5khz and f_high 
ind_freq = (fc<=f_high) & (fc>=f_low); % relevent freqs
% ind_freq = [ind_freq(1)-1, ind_freq, ind_freq(end)+1];
fc = fc(ind_freq); % pick freqs between f_low,f_high
a = 1; % atten
n_order = 8; % order of the filter (4)
t = [0:ts:response_time]; % time vector
f_in = [f_low: fft_resolution: f_high];
rise_ratio = 0.8; % the percentage of the rise in amplituude to measure
coeff_fix = 10*(250000/fs); % measured to improve rise time esimation 

test_rise_flg =0;

% init varaibles
impulse_resp = zeros(numel(fc), numel(t));
delay_times = zeros(1,numel(fc));
P1 = zeros(numel(fc), floor(numel(t)/2)+1);
fft_resp = zeros(numel(fc), numel(f_in));
cross_matrix = zeros(numel(fc));

plot_flag = 0;
if plot_flag
    figure()
    ax_imp = subplot(2,1,1);
    xlabel(ax_imp, 'time');
    hold(ax_imp,'off')
    
    ax_fft = subplot(2,1,2);
    xlabel(ax_fft, 'freq');
    ylabel(ax_fft, 'response');
    hold(ax_fft,'on')
    
    set(gcf, 'Units','inches', 'Position', [ 2 1 9 6], 'PaperPositionMode','auto');
    ax_fft.Position =  [0.13   0.08    0.7750    0.4];
    ax_imp.Position =  [0.13   0.56    0.7750    0.4];
end % plot_flag


%% the filter bank in each fc
% rect_signal = ones(1,floor(numel(t)/5));

for k =1:numel(fc)
    b = 0.15*fc(k); % bandwidth of the filter, relative to fc
    % the gamma filter
    gt = a*t.^(n_order-1) .* exp(-2*pi*b*t) .* cos(2*pi*fc(k)*t);
    % the impulse response
%     lpf_gt = envelope(gt);
    impulse_resp(k,:) = gt ./ max(abs(gt));% lpf_gt/max(abs(lpf_gt));
%     [~, delay_times(k)] =  max(impulse_resp(k,:));
  


    %% XXX clculate rise-times of the Pure Impuse respose for each filter

        % delay time to maximum
        [pks, ipks] = findpeaks(envelope(impulse_resp(k,:)), 'MinPeakDistance',20,  'MinPeakHeight',0.5, 'NPeaks', 1);
        delay_times(k) = ts*ipks(1);
            
    %% end XXX
    %% fft_response - It's not relevant for acoustic but still 
    resp_fft = fft(gt); %fft(impulse_resp(k,:));
    L = numel(resp_fft);
    resp_fft = abs(resp_fft)/L;
    P1(k,:) = resp_fft(1:floor(L/2)+1);
    P1(k, 2:end-1) = 2*P1(k,2:end-1);
    P1(k,:) = P1(k,:)/max(P1(k,:));
    f = fs*(0:(L/2))/L;
    fft_resp(k,:)= interp1(f,P1(k,:),f_in);
    cross_matrix(k,:) = interp1(f,P1(k,:),fc); % the matrix of the cross talk between fcs
    
    if plot_flag && mod(k,5) == 1  
        plot(ax_imp, t, impulse_resp(k,:),'DisplayName', [num2str(round(fc(k))), ': ', num2str(t(ipks)*1e3,2), 'ms' ])
%         plot(ax_imp, t(ipks), pks, '*r')
%         text(ax_imp, t(ipks), pks*1.05, num2str(round(fc(k))) );
        plot(ax_fft, f, P1(k,:),'DisplayName', num2str(fc(k),2));
%         plot(ax_fft, f_in, fft_resp(k,:),'o')
        legend (ax_imp)
    end % plot_flag
end % for k


AllFiltersParams = struct( ...
    'filter_fc', fc, ...
    'filter_impulse_resp', impulse_resp, ...
    'filter_delay_times', delay_times, ...
    'fft_response', fft_resp, ...
    'filter_timevec', t, ... 
    'fft_freqs', f_in, ...
    'cross_matrix', cross_matrix);
