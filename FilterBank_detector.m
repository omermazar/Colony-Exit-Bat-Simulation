function [fc_detect_struct, integrated_detect] = FilterBank_detector(rx_strct, FilterBankRxMat, JammingMat, FilterBank, AllParams, det_method)
% inputs:
% rx_strct = BAT(BatNumRx).FilterBank_Pulses_strct(PulseNum); % struct of the rx signals in this calls
% general parameters
% FilterBankRxMat: the Rx matrix of setection in us
% det_method : 'Fcs'- histogram of detections in channels , 'Integrated' - ingration of the channels and than decision 

%% general
LinTH = 10.^(AllParams.BatSonarParams.PulseDetectionTH/10); % linear th
Ts = 1/AllParams.BatSonarParams.FilterBank_Fs; 
fs_ts = AllParams.SimParams.SampleTime * AllParams.BatSonarParams.FilterBank_Fs;
NoiseLevel = 10.^(AllParams.BatSonarParams.NoiseLeveldB/10);
plot_fc_flag = 0;

%% parameters for findpeaks
% calculate the duration of each transmitted fc
fc_start = [0, find(diff( rx_strct.TxActiveFcs_idx))']; 
fc_dur = diff([fc_start, numel(rx_strct.TxActiveFcs_idx)]);
mean_dur = mean(fc_dur);

min_peak_sense = 0.5*NoiseLevel; % for the peak of the rise time (srivative) To test 2 is tooo high 
                    % 0.5 - good detection , no interference
min_peak_ct = 8*NoiseLevel; % 5 % the th for peak detection  
rise_time_th = 0.9; % the rise time is (rise_time_th*peak)

min_diff_peaks = max(20, mean_dur/3);  % 50; % samples (200 micro) % 20 XXXX
grad_param = 10; % gain parameter for gradient
CrossTalkFlag = FilterBank.CrossTalkFlg;
short_pulses_correction = 10; % samples

% max_rise_samples= round(1e-3/Ts); % 1ms
% rise_slope_ratio = 0.2;

%% Input data
% the relevant times (in samples (us)) to detect signals
% % tt = (rx_strct.TxStartTime) : (rx_strct.TxStartTime + rx_strct.PulseIPI); % XXX
tt = 0: rx_strct.PulseIPI;
tt_ind = 1:numel(tt);
% the relevant fcs ('filter-bank-channels')
unique_active_fcs = unique(rx_strct.TxActiveFcs,'stable');
unique_active_fcs_idx = unique(rx_strct.TxActiveFcs_idx,'stable');
max_fcs_samples = max([125, max(diff(rx_strct.Tx_Fcs_StartTimes))]);


%% init output for decision
integrated_detect = struct(...
    'int_signal', [], ...
    'pks', [], ...
    'delays', [], ...
    'wdths',[]);

%% THE DETECTION %%%%%%%%%%%%%%%%%%%%%%%
switch det_method
    
    case 'Integrated'
        %% Integerated detection
%         int_ir_delay = 0;
        % the error is biased as function of the pulse -duration
        m = 0.04;
        b = 30;
        int_ir_delay = m * rx_strct.PulseDuration +b; 
%         int_ir_delay = mean(FilterBank.cross_talk_delays(unique_active_fcs_idx));
        noise_mat = 1/sqrt(2)*NoiseLevel*randn(numel(unique_active_fcs_idx), size(FilterBankRxMat,2));
        mat_sig_in = FilterBankRxMat(unique_active_fcs_idx,:) + ...
            JammingMat(unique_active_fcs_idx,:) ...
            +  noise_mat;
        integrated_pp = fb_integrate_filters(rx_strct, mat_sig_in);
        
        [pks_int, locs_int] = findpeaks( integrated_pp, 'MinPeakDistance',min_diff_peaks  , ...
            'MinPeakHeight', 5*NoiseLevel, 'MinPeakProminence', min_peak_sense); % [pks_int, locs_int, wdths_int, proms]
        
        k_factor = numel(unique_active_fcs);
        
        %%%% XXXXX %%%%%
        %%% while emitting a pulsethe bat cannot detect
        
        tx_idx = locs_int < 0.6*rx_strct.PulseDuration;
        if sum(tx_idx) > 0
            locs_int = locs_int(~tx_idx);
            pks_int = pks_int(~tx_idx);
        end %if sum(not_tx_idx) > 0
        %%%% XXXXX %%%%%
        
        % debug
        integrated_detect.int_signal = integrated_pp; % for debug and plot mainly
        integrated_detect.pks = pks_int;
        integrated_detect.delays = tt(locs_int); % + rx_strct.TxStartTime;
        integrated_detect.wdths = []; % wdths_int;
        integrated_detect.proms = []; %proms;
        
        % outputs
        fc_detect_struct.active_fc = unique_active_fcs; 
        fc_detect_struct.active_fc_idx = unique_active_fcs_idx; 
        fc_detect_struct.tx_fcs_start_time = []; 
        fc_detect_struct.num_of_detections = numel(pks_int); 
        fc_detect_struct.raw_detections = locs_int; 
        fc_detect_struct.detection_delays = fc_detect_struct.raw_detections -int_ir_delay; 
        ind_neg_idx = find(fc_detect_struct.detection_delays <=0);
        if ~isempty(ind_neg_idx)
             fc_detect_struct.detection_delays(ind_neg_idx) =1;
        end %  ~isempty(ind_neg_idx)
        fc_detect_struct.rise_time_det = fc_detect_struct.detection_delays + rx_strct.TxStartTime;
        fc_detect_struct.rx_power_db  = 10*log10(pks_int/k_factor);
        
        
    case 'Fcs'
        %% post process fcs
        noise_vec = 1/sqrt(2)*NoiseLevel*randn(1,size(FilterBankRxMat,2));
        % init output for decision
        fc_detect_struct(numel(unique_active_fcs_idx))= struct(...
            'active_fc', [],...
            'active_fc_idx', [], ...
            'tx_fcs_start_time',[], ...
            'num_of_detections', [], ...
            'rise_time_det', [], ...
            'detection_delays', [], ...
            'raw_detections', [], ...
            'rx_power_db',[]);

        post_process_flag= 1;
        if post_process_flag
            pp_method = 4;
            pp_gain = 4; % adding 6 db for common ground with not filter bank receivers
            [FilterBankRxMat_pp] = fb_postprocess(rx_strct, FilterBankRxMat, FilterBank, pp_method, pp_gain);
            
            % integrated estimation of the delay (equivalnt to semi correlation
            % %     integrated_pp = smooth(integrated_pp + noise_vec, numel(unique_active_fcs_idx)); %#ok<NASGU>
%             [pks_int, locs_int, wdths_int] = findpeaks( integrated_pp, 'MinPeakDistance',min_diff_peaks  , ...
%                 'MinPeakHeight', 5*NoiseLevel, 'MinPeakProminence', min_peak_ct);
%             % %
%             integrated_detect.int_signal = integrated_pp; % for debug and plot mainly
%             integrated_detect.pks = pks_int;
%             integrated_detect.delays = tt(locs_int); % + rx_strct.TxStartTime;
%             integrated_detect.wdths =  wdths_int;
            % %
        end % if post_process_flag
        
        %% detection in each fc
        
        kk = 1;
        for k_idx = unique_active_fcs_idx'
            fc_detect_struct(kk).active_fc = unique_active_fcs(kk);
            fc_detect_struct(kk).active_fc_idx = k_idx;
            % thestart-time of the transmitted signal
            tx_fcs_start_time = rx_strct.TxTimes_us(find(rx_strct.TxActiveFcs_idx == k_idx, 1));
            fc_detect_struct(kk).tx_fcs_start_time = tx_fcs_start_time;
            % the rx signal in the fc-channel:
            %%% ________________%%%
            if post_process_flag
                rx_vect = FilterBankRxMat_pp(kk,tt_ind) + JammingMat(k_idx,tt_ind); % check
            else % post_process_flag
                rx_vect = FilterBankRxMat(k_idx,tt_ind) + JammingMat(k_idx,tt_ind); % check
            end % if post_process_flag
            %     rx_vect = FilterBankRxMat(k_idx,tt_ind) + JammingMat(k_idx,tt_ind); % check
            %%% ________________%%%
            
            % add noise to the signal, signal below noise is not relevealt
            
            rx_vec = rx_vect + noise_vec; %% max(rx_vect, NoiseLevel*ones(size(rx_vect)) );
            
            max_rx = max(rx_vec);
            
            % condition #1: detection threshold in each fc
            is_detected_th = max_rx > LinTH;
            if is_detected_th
                %%% find peaks of the response
                if CrossTalkFlag == 0
                    dy =  gradient(rx_vec, 1/grad_param); % find rise times
                    [dypks,idt] = findpeaks(dy, 'MinPeakDistance',min_diff_peaks, ...
                        'MinPeakHeight', NoiseLevel, 'MinPeakProminence', min_peak_sense);
                    % only positive derivative peaks
                    pos_pks = dypks>0;
                    dypks = dypks(pos_pks);
                    idt = idt(pos_pks);
                    fc_ir_delay = FilterBank.filter_delay_times(k_idx);
                    
                elseif CrossTalkFlag == 1
                    %              [dypks,idt, wdths] = findpeaks(rx_vec, 'MinPeakDistance',min_diff_peaks, ...
                    %                   'MinPeakHeight', 5*NoiseLevel, 'MinPeakProminence', min_peak_ct);
                    [dypks,idt, ~] = findpeaks(rx_vec, 'MinPeakDistance',min_diff_peaks, ...
                        'MinPeakHeight', 5*NoiseLevel, 'MinPeakProminence', min_peak_ct);
                    fc_ir_delay = FilterBank.cross_talk_delays(k_idx);
                    
                end % if CrossTalkFlag == 0
                fc_detect_struct(kk).num_of_detections = numel(dypks);
                
                % condition #2: detection in each fc
                if fc_detect_struct(kk).num_of_detections > 0
                    
                    %% estimate rise time from peaks
                    % %             delay_coeff =  0 ; % FilterBank.rise_time_coeff(k_idx);
                    % %             rise_times = round(idt + delay_coeff); % samples
                    rise_times= zeros(size(dypks));
                    
                    if CrossTalkFlag
                        for m = 1:numel(dypks)
                            %                     il = max(1, idt(m)-round(wdths(m)) );
                            %                     ir = find(rx_vec(il:idt(m)) >= rise_time_th*dypks(m), 1);
                            il = idt(m);
                            if il > fc_dur(kk)               %  wdths(m)
                                min_ind = round(il-fc_dur(kk));   % round(il-wdths(m));
                            else
                                min_ind = 1;
                            end % if %il > width(m)
                            ir = find(rx_vec(il:-1:min_ind ) <= rise_time_th*dypks(m), 1);
                            if ~isempty(ir)
                                rise_times(m) = il-ir+1;
                            else % if ~isempty(ir)
                                %                         rise_times(m) = max([1, round(il-wdths(m)/2) ]);   % protect start of signal
                                rise_times(m) = max([1, round(il-fc_dur(kk)/2) ]);   % protect start of signal
                            end % if ~isempty(ir)
                        end % for m
                        
                    end % if CrrossTalkFlag
                    
                    
                    %%
                    
                    %% output of the channel
                    %             delays_detected = tt(rise_times) + rx_strct.TxStartTime- ...
                    %                 ( tx_fcs_start_time + FilterBank.filter_delay_times(k_idx) );
                    delays_raw = tt(rise_times  ) + rx_strct.TxStartTime- ...
                        ( tx_fcs_start_time );
                    delays_detected = delays_raw- fc_ir_delay;%  delays_raw - fc_ir_delay;
                    
                    %%% correction for final buzz short pulses
                    if rx_strct.PulseDuration <= fs_ts
                        delays_detected = delays_detected + short_pulses_correction;
                    end % if rx_strct.PulseDuration <= fs_ts
                    % output
                    % %             fc_detect_struct(kk).rise_time_det = tt(rise_times); % XXX tt+ txtime(1
                    fc_detect_struct(kk).rise_time_det = tt(rise_times) + rx_strct.TxStartTime;
                    fc_detect_struct(kk).detection_delays = delays_detected;
                    fc_detect_struct(kk).raw_detections = delays_raw;
                    
                    
                    % calculate rx_power of the pulses relwvant for risetimes
                    if numel(rise_times) > 1
                        diff_rt = [diff(rise_times), max(diff(rise_times))];
                    else % if numel(rise_times) > 1
                        diff_rt = max_fcs_samples;
                    end % if numel(rise_times) > 1
                    rx_power = zeros(size(rise_times));
                    for k_pk = 1:numel(rise_times)
                        try
                            max_ind = min(numel(tt), rise_times(k_pk)+diff_rt(k_pk) );
                            if rise_times(k_pk) < max_ind
                                rx_power(k_pk)= max(rx_vec(rise_times(k_pk): max_ind) );
                            end % if rise_timr
                        catch
                            kk
                            k_pk
                        end
                    end
                    fc_detect_struct(kk).rx_power_db = 10*log10(abs(rx_power));
                    
                    
                end % if fc_detect_stuct(kk).detection_delays > 0
            else % not detected in this channel
                fc_detect_struct(kk).num_of_detections = 0;
            end % if is_detected_th
            if plot_fc_flag
                plot_detector_analyze
            end %if plot_vc_flag
            kk= kk+1;
        end % for k_idx
        
end % switch det_method
