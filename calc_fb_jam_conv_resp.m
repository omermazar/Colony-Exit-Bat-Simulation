function [rx_mat_updated, t_out]= calc_fb_jam_conv_resp( jam_strcts, rx_mat, rx_fb_strct, ...
    start_time, end_time, FilterBank, fs_ts, sample_limit)

% [rx_mat_updated]= calc_fb_conv_resp( rx_strct, rx_mat, FilterBank, calc_flag)
% The function calcultes the conolution beween the received signals and the
% impulse response of the filter bank 
%
% output: 
%   rx_mat: Matrix of the intensities in each 'channel' and times
%   rx_vec: vector of the detection in times (1 if it is a signal, 0 if
%   not)
% inputs:
%    jam_strcts: a sruct of the relevant recpetion from other bats, in jam_strcts(Times):
%   struct with fields:
%             Times: 152
%     Times_precise: 152.4374
%             Freqs: 49.7090
%             Power: 9.0996e+03
%          BATTxNum: 2
%     BatTxPulseNum: 1

%   rx_mat: the matrix to be updated, the rows in rx_mat are the rx signal in each fft signal before convolution
%   FilterBank: a struct with the impulse responses of the FilterBank       
%   struct with fields:
%               filter_fc: [1×63 double]
%     filter_impulse_resp: [63×126 double]
%      filter_delay_times: [1×63 double]
%            fft_response: [63×71 double]
%          filter_timevec: [1×126 double]
%               fft_freqs: [1×71 double]
%   fs_ts =  Fs * SampleTimet : he ratio between the simulation sample_time and the hih
%       reolustion filter-bank's sample-time
%   calc_flag: 'Rx' for own echoes, 'Tx' for intereference

CrossTalk_flag = FilterBank.CrossTalkFlg;
max_t_out = sample_limit; %end_time - start_time +1;

max_time_rx = max_t_out; % min(jam_strct.PulseIPI, max_t_out); % max index to avoid ambiguiuty

orig_jam_BatTx = vertcat(jam_strcts.BATTxNum);
    
tt= start_time: end_time;
tt_ind = 1:numel(tt);
conv_time = numel(tt) + numel(FilterBank.filter_impulse_resp(1,:))-1;
t_out = tt(1): min(max_t_out, tt(1) + conv_time-1);

%% calculate impulse response for each row in rx_mat 
rx_mat_updated = zeros(numel(FilterBank.filter_fc), numel(t_out) ); 
if CrossTalk_flag
    % sum of all active freq with cross talk
    cur_freq_resp = FilterBank.fft_response * rx_mat;
    % calculate the response for relevant fcs only 
    rx_fcs_idx = unique(rx_fb_strct.TxActiveFcs_idx);
    for k_row= rx_fcs_idx'  % 1:size(cur_freq_resp,1)
        ir = conv(cur_freq_resp(k_row,:), FilterBank.filter_impulse_resp(k_row,:))/(2*sqrt(fs_ts));
        if numel(ir) ~= size(rx_mat_updated,2)
            ir = ir(1:size(rx_mat_updated,2));
        end % if numel
        rx_mat_updated(k_row,:) = ir;
    end
    
else %  if CrossTalk_flag - no CrossTalk
    for k_freq_in= 1:size(rx_mat,1)
        is_empty = sum(rx_mat(k_freq_in,:)) == 0;
        if is_empty % chech if ther is any signal in that channel
            continue
        else % if is_empty
            % The rx signal
            rx_ch = rx_mat(k_freq_in,:);
            % the relvant filter
            v_fcs = 1:numel(FilterBank.filter_fc);
            fcs_idx =  round(interp1(FilterBank.filter_fc, v_fcs, FilterBank.fft_freqs(k_freq_in)) ); % the index of the fcs response
            if any(isnan(fcs_idx))
                fcs_idx(isnan(fcs_idx)) = numel(FilterBank.filter_fc); % the maximum freq
                fft_idx(isnan(fft_idx)) = numel(FilterBank.fft_freqs);
            end % if isnan
            cur_freq_resp = FilterBank.fft_response(fcs_idx,k_freq_in); % vector with the freq response in the active_fcs in time
            % the imules response
            ir = conv(cur_freq_resp.*rx_ch, FilterBank.filter_impulse_resp(fcs_idx,:))/(2*sqrt(fs_ts)); % XXX
            ir= ir(1:numel(t_out)); % end of time
            rx_mat_updated(fcs_idx,:) = ir; % the response is in the relevant fc channel
        end %if is_empty
        
    end % k_freq_in= 1:size(rx_mat,1)
end %  if CrossTalk_flag

        

        %%
        %% figure
        plot_flag = 0;
        if plot_flag
            figure
            title(['echo', num2str(k_echo)])
            xlabel('time')
            ylabel('channel ir')
            hold on
            plot(tt, active_fcs/1000)
            plot(tt, 10*log10(rx_ch_tot))
            
            for k_idx = 1:numel(rx_fcs_idx)
                k_fcs = rx_fcs_idx(k_idx);
                plot(t_out, 10*log10(echo_mat_in(k_fcs,:)),'LineWidth',1.5);
                text( min(tt)+100, -20-k_idx*5, num2str(FilterBank.filter_fc(k_fcs)), 'FontSize', 8)
            end % for k_idx
            
        end % if plot_flag
        
    end % if rx_strct.RxEchosFromPrey(k_echo).Rx_times < max_time_rx

  % the         



    