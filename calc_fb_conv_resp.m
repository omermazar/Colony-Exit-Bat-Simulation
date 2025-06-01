function [rx_mat_updated, rx_track_updated]= calc_fb_conv_resp( rx_strct, rx_mat, rx_track_in, ...
    start_time, end_time, FilterBank, fs_ts)

% [rx_mat_updated]= calc_fb_conv_resp( rx_strct, rx_mat, FilterBank, calc_flag)
% The function calcultes the conolution beween the received signals and the
% impulse response of the filter bank 
%
% output: 
%   rx_mat: Matrix of the intensities in each 'channel' and times
%   rx_vec: vector of the detection in times (1 if it is a signal, 0 if
%   not)
% inputs:
%    rx_strct: a sruct of the relevant recpetion:
%   struct with fields:
%              PulseNum: 1
%            PulsePower: 1.0000e+11
%            TxTimes_us: [1×1750 double]
%           TxStartTime: 6500
%         PulseDuration: 1750
%          TxFrequncies: [1×1750 double]
%           TxActiveFcs: [1×1750 double]
%       TxActiveFcs_idx: [1750×1 double]
%     Tx_Fcs_StartTimes: [6500 6625 6750 6875 7000 7125 7250 7375 7500 7625 7750 7875 8000 8125]
%       RxEchosFromPrey: [1×10 struct]

%   Rx_mat: the matrix to be updated
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

CrossTalk_flag = FilterBank.CrossTalkFlg  ;
plot_flag = 0;
if plot_flag
    figure ()
end % if plot_flag

tx_update_flag = 0;
% max_t_out = size(rx_mat,2); % the end of the relevant duration to calc simulation XXX
max_t_out = end_time - start_time +1;
ref_star_time = start_time;

rx_mat_updated = rx_mat;
rx_track_updated = rx_track_in;
% max_time_rx = min(rx_strct.TxStartTime + rx_strct.PulseIPI, max_t_out); % max index to avoid ambiguiuty
max_time_rx = min(rx_strct.PulseIPI, max_t_out); % max index to avoid ambiguiuty


    
%% update rx_mat with Tx Values
if tx_update_flag
    idx_mat = sub2ind(size(rx_mat_updated), rx_strct.TxActiveFcs_idx', rx_strct.TxTimes_us+1 - start_time);
    rx_mat_updated(idx_mat) = rx_mat_updated(idx_mat) + rx_strct.PulsePower;
end % if tx_update_flag


%% calculate impulse response in each fc ('channel') for each echo in pulse
for k_echo = 1:length(rx_strct.RxEchosFromPrey)
    
    %  avoid ambiguity
    if rx_strct.RxEchosFromPrey(k_echo).Rx_times(1) - start_time < max_time_rx-1
        tt = rx_strct.RxEchosFromPrey(k_echo).Rx_times - start_time + 1 ; % vector of the recived times of the signal
        tt = round( tt(1) : min(tt(end), max_time_rx) );
        tt_ind = 1:numel(tt);
        % the 'sampled' frequencies of the signal
        v_fft = 1:numel(FilterBank.fft_freqs);
        fft_idx = round(interp1(FilterBank.fft_freqs, v_fft, rx_strct.RxEchosFromPrey(k_echo).Rx_freqs(tt_ind)))';
%         [fft_idx,~] = dsearchn(FilterBank.fft_freqs', rx_strct.RxEchosFromPrey(k_echo).Rx_freqs(tt_ind)');
%         un_fft_idx = unique(fft_idx, 'stable');
%         active_fft_freqs = FilterBank.fft_freqs(fft_idx);
        
% the matrix of the response in each channel
        cur_freq_resp = zeros(numel(FilterBank.filter_fc), numel(tt));
        v_fcs = 1:numel(FilterBank.filter_fc);
        fcs_idx =  round(interp1(FilterBank.filter_fc, v_fcs, rx_strct.RxEchosFromPrey(k_echo).Rx_freqs(tt_ind)))';
        if any(isnan(fcs_idx))
            fcs_idx(isnan(fcs_idx)) = numel(FilterBank.filter_fc); % the maximum freq
            fft_idx(isnan(fft_idx)) = numel(FilterBank.fft_freqs); % the maximum freq
        end % if isnan
%         [fcs_idx,~] = dsearchn(FilterBank.filter_fc', rx_strct.RxEchosFromPrey(k_echo).Rx_freqs(tt_ind)');
        active_fcs = FilterBank.filter_fc(fcs_idx);
        un_fcs_idx = unique(fcs_idx, 'stable');
        
        %%% for  the conveloution of each channel with each impulse response
        % each raw is an fc-channel
        conv_time = numel(tt) + numel(FilterBank.filter_impulse_resp(1,:))-1;
        echo_mat_in = zeros(numel(FilterBank.filter_fc), conv_time  );
        t_out = tt(1): min(max_t_out, tt(1) + conv_time-1);
        
        
        if CrossTalk_flag
            %% XXX 

            rx_ch_tot = zeros([numel(FilterBank.fft_freqs), length(tt)]); % the rx signal power in time
            
            mat_ct_idx = sub2ind( size(rx_ch_tot), fft_idx, tt_ind');
% %             mat_ct_idx = sub2ind(size(cur_freq_resp), fft_idx, [1:length(tt)]');
            rx_ch_tot(mat_ct_idx) = rx_strct.RxEchosFromPrey(k_echo).Rx_powers_lin(tt_ind);
            cur_freq_resp = FilterBank.fft_response * rx_ch_tot;
            % conv in all_ffts
            for k_raw= un_fcs_idx'  % 1:size(cur_freq_resp,1)
                ir = conv(cur_freq_resp(k_raw,:), FilterBank.filter_impulse_resp(k_raw,:));
%                 ir = conv(cur_freq_resp(k_raw,:).*rx_ch_tot(k_raw,:), FilterBank.filter_impulse_resp(k_raw,:));
                echo_mat_in(k_raw,:) = echo_mat_in(k_raw,:)+ ir/(2*sqrt(fs_ts));
            end
            %%
        else % if CrossTalk_flag - no CrossTalk
            fcs_fft_idx = sub2ind(size(FilterBank.fft_response), fcs_idx,  fft_idx);
            cur_freq_resp_vec = FilterBank.fft_response(fcs_fft_idx); % vector with the freq response in the active_fcs in time
            mat_idx = sub2ind(size(cur_freq_resp), fcs_idx, [1:length(tt)]');
            cur_freq_resp(mat_idx) = cur_freq_resp_vec;
            rx_ch_tot = zeros(1, size(cur_freq_resp,2)); % the rx signal power in time
            
            % conv in active_fcs
            for k_ch = 1:numel(un_fcs_idx)
                k_fcs = un_fcs_idx(k_ch);
                fc_time_idx = fcs_idx == k_fcs;
                rx_ch = zeros(1, size(cur_freq_resp,2));
                
                % the ampilude as function of time in the channel (no croesstoalk
                rx_ch(fc_time_idx) = rx_strct.RxEchosFromPrey(k_echo).Rx_powers_lin(fc_time_idx);
                rx_ch_tot = rx_ch_tot+ rx_ch;
                
                ir = conv(cur_freq_resp(k_fcs,:).*rx_ch, FilterBank.filter_impulse_resp(k_fcs,:));
                echo_mat_in(k_fcs,:) = ir/(2*sqrt(fs_ts));
            end
        end % if if CrossTalk_flag
        
        
        
        
        %% update rx_mat
        if t_out(end) == max_t_out
            echo_mat_in = echo_mat_in(:,1:numel(t_out));
        end %         if t_out(end) == max_t_out

        rx_mat_updated(:,round(t_out)) = rx_mat_updated(:,round(t_out)) + echo_mat_in;
        rx_track_updated(rx_strct.RxEchosFromPrey(k_echo).target_idx,   round(t_out)) = rx_strct.RxEchosFromPrey(k_echo).target_idx;
        %%
        %% figure

        if plot_flag
%             figure
            hold on
            title(['echo', num2str(k_echo)])
            xlabel('time')
            ylabel('channel ir')
            hold on
            plot(tt, active_fcs/1000,'m.')
            plot(tt, 10*log10(rx_ch_tot),'b','LineWidth',1.5)
            
%             for k_idx = 1:numel(un_fcs_idx)
%                 k_fcs = un_fcs_idx(k_idx);
%                 plot(t_out, 10*log10(echo_mat_in(k_fcs,:)),'LineWidth',1.5);
%                 text( min(tt)+100, -20-k_idx*5, num2str(FilterBank.filter_fc(k_fcs)), 'FontSize', 8)
%             end % for k_idx
            
        end % if plot_flag
        
    end % if rx_strct.RxEchosFromPrey(k_echo).Rx_times < max_time_rx
end % for k_echo

tt_out = [1:(size(rx_mat_updated,2))];

if plot_flag
    
    for k_idx = 1:numel(un_fcs_idx)
        k_fcs = un_fcs_idx(k_idx);
        plot(tt_out, 10*log10(abs(rx_mat_updated(k_fcs,:))));
        text( min(tt)+100, -20-k_idx*5, num2str(FilterBank.filter_fc(k_fcs)), 'FontSize', 8)
    end % for k_idx
    
end % if plot_flag


    