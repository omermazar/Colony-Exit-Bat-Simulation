function [AllFiltersParams]  = ...
    FilterBank_Resp(num_of_filters, fs, f_low, f_high, response_time, fft_resolution)

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
%   f_low = lowest freq in thc filterbank (center freq) (default: 10e3)
%   f_high - highest freq in the filter bank (default(80e3)
%   response_time - the maximum time of the impluse response (500e-6)
%   fft_resolution - the resolution of the fft output (1e3)
% Default Params:

%% general parameters
% num_of_filters = 80;
% fs = 0.1e6;
% f_low = 32e3;
% f_high = 80e3;
ts= 1/fs; 
k_filt = 1:num_of_filters;
k_norm = num_of_filters/ log2(f_high/5703);
fc = 5702*2.^(k_filt/k_norm); % total of num_of_filters between ~5khz and f_high 
ind_freq = (fc<=f_high) & (fc>=f_low); 
% ind_freq = [ind_freq(1)-1, ind_freq, ind_freq(end)+1];
fc = fc(ind_freq); % pick freqs between f_low,f_high
a = 1; % atten
n_order = 4; % order of the filter
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
    ax_imp = gca;
    xlabel(ax_imp, 'time');
    hold(ax_imp,'on')
    
%     figure()
%     ax_fft = gca;
%     xlabel(ax_fft, 'freq');
%     ylabel(ax_fft, 'response');
%     hold(ax_fft,'on')
end % plot_flag
%%
%% input signal - Impulse 
% t_in = t; % [0:sample_time:10e-3];
% sig_in = zeros(size(t_in));
% sig_in(1) = 1;

%% the filter bank in each fc
cross_talk_delays= zeros(size(fc)); %%%
rect_signal = ones(1,floor(numel(t)/5));

for k =1:numel(fc)
    b = 0.15*fc(k); % bandwidth of the filter, relative to fc
    % the gamma filter
    gt = a*t.^(n_order-1) .* exp(-2*pi*b*t) .* cos(2*pi*fc(k)*t);
    % the impulse response
    lpf_gt = envelope(gt);
    impulse_resp(k,:) = lpf_gt/max(abs(lpf_gt));
%     [~, delay_times(k)] =  max(impulse_resp(k,:));
   
    
    %% XXX clculate rise-times
    
% %         [max_resp, imax] =  max(impulse_resp(k,:));
% %         rise_value= rise_ratio*max_resp;
% %         [~,delay_times(k)] = min( abs(impulse_resp(k,1:imax) - rise_value));
% %         %     delay_times(k) = min(iv);
        ref_sig = conv(impulse_resp(k,:),rect_signal)/numel(rect_signal); % rect_signal)/25;
        dy = gradient(ref_sig,0.1);
%         dy = gradient(impulse_resp(k,:),0.1);         [dypks, idy] = findpeaks(dy, 'MinPeakDistance',20,  'MinPeakHeight',0.2);
        [~, idy] = findpeaks(dy, 'MinPeakDistance',20,  'MinPeakHeight',0.2);
        delay_times(k) = idy(1);
        % delay time to maximum
        [~, ipks] = findpeaks(ref_sig, 'MinPeakDistance',20,  'MinPeakHeight',0.5);
        cross_talk_delays(k) = ipks;  % max(-idy+1, delay_times(k)- idy- coeff_fix);
            
    %% end XXX
    %% fft_response
    resp_fft = fft(gt); %fft(impulse_resp(k,:));
    L = numel(resp_fft);
    resp_fft = abs(resp_fft)/L;
    P1(k,:) = resp_fft(1:floor(L/2)+1);
    P1(k, 2:end-1) = 2*P1(k,2:end-1);
    P1(k,:) = P1(k,:)/max(P1(k,:));
    f = fs*(0:(L/2))/L;
    fft_resp(k,:)= interp1(f,P1(k,:),f_in);
    cross_matrix(k,:) = interp1(f,P1(k,:),fc); % the matrix of the cross talk between fcs
    
    if plot_flag  %% & ~mod(k,4)
        
        plot(ax_imp, t, impulse_resp(k,:))
        plot(ax_imp, t(), impulse_resp(k,delay_times(k)), 'o')
% %         plot(ax_imp, t, ref_sig(1:length(t)), 'r')
        plot(ax_imp, t(cross_talk_delays(k)), ref_sig(cross_talk_delays(k)), '*r')
% %         figure
% %         hold on
% %         title(['fcs= ', num2str(fc(k))])
% %         plot(impulse_resp(k,:),'b')
% %         plot(delay_times(k), impulse_resp(k,delay_times(k)), 'ob')
% %         plot( ref_sig(1:length(t)), 'r')
% %         plot(cross_talk_delays(k), ref_sig(cross_talk_delays(k)), '*r')
%             plot(ax_imp, t(ipks), pks, 'dr')

%         plot(ax_fft, f, P1(k,:));
%         plot(ax_fft, f_in, fft_resp(k,:),'o')
        
    end % plot_flag
end % for k

AllFiltersParams = struct( ...
    'filter_fc', fc, ...
    'filter_impulse_resp', impulse_resp, ...
    'filter_delay_times', delay_times, ...
    'fft_response', fft_resp, ...
    'filter_timevec', t, ... 
    'fft_freqs', f_in, ...
    'cross_talk_delays', cross_talk_delays, ...
    'cross_matrix', cross_matrix);
