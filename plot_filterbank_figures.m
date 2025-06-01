% function [] = plot_filterbank_figures(AllParams,CurrEchosFromPreyStruct, rx_strct, tt, rx_vec)

figure()
hold on
amp_ax = gca;
% title(['ampltude figures, fc= ', num2str(round(fc_detect_struct(kk).active_fc))])
%% rx-power plot
%% input signal
allpower_flg =1;
if allpower_flg
%     tx_power = 10.^(BAT(BatNum).TransmittedPulsesStruct(BAT(BatNum).CurrPulseNum).PulsePower/10);
    tx_power = 10.^(BAT(BatNum).TransmittedPulsesStruct(CurrentPulseNum).PulsePower/10);
    for k=1:AllParams.SimParams.TotalPreysNumber
        
%         plot(amp_ax, CurrEchosFromPreyStruct.EchoDetailed(k).EchosnTimesFull*125, ...
%             abs(CurrEchosFromPreyStruct.EchoDetailed(k).EchoAttenuationFull*tx_power) ,'o-')
        plot(amp_ax, CurrEchosFromPreyStruct.EchoDetailed(k).EchosnTimesFull*125, ...
            10*log10(abs(CurrEchosFromPreyStruct.EchoDetailed(k).EchoAttenuationFull*tx_power)) ,'o-')
%         %         plot(amp_ax, rx_strct.RxEchosFromPrey(k).Rx_times, ...
        %             10*log10(rx_strct.RxEchosFromPrey(k).Rx_powers_lin), '.')
        
    end %for k
    plot(amp_ax, amp_ax.XLim, zeros(size(amp_ax.XLim)), 'r--') % reference TH
end % if allpower

 for k=1:AllParams.SimParams.TotalPreysNumber
%         plot(CurrEchosFromPreyStruct.EchoDetailed(k).EchosnTimesFull*125,...
%             CurrEchosFromPreyStruct.EchoDetailed(k).EchosFreqsFull ,'o-')
        plot(rx_strct.RxEchosFromPrey(k).Rx_times, round(rx_strct.RxEchosFromPrey(k).Rx_freqs/1000), '.')
 end

 %%  input ro Filterbank detector
 figure
 hold on
 title('signals')
 
 tt_start_time = rx_strct.TxStartTime;

 rx = FilterBankRxMat(unique_active_fcs_idx,tt_ind); 
    %  rx = FilterBankRxMat(k_idx:,tt_ind);
 jam = JammingMat(unique_active_fcs_idx,tt_ind);
 % jam = JammingMat(k_idx:,tt_ind);
%  all channels
 plot( tt+tt_start_time,10*log10(rx+0.01),'k')
 plot( tt+tt_start_time,10*log10(jam+0.01),'r')
%  plot( tt+tt_start_time,10*log10(rx_vec),'.')
 
 %% output total
 plot([fc_detect_struct.rise_time_det], [fc_detect_struct.rx_power_db],'*r')
 
 rx_vect2= rx;
 int_lvl = max(jam, NoiseLevel*ones(size(rx_vect)));
 plot( tt+tt_start_time,10*log10(int_lvl),'or')
 rx_vec2 = max(rx, int_lvl);
 
 %% integrated output
 figure
 hold on
 title(['integrated detetor and echoes, pulse: '])
 
 %%%% if ruuning from main
 main_run = 0;
 if main_run
     rx_strct = BAT(BatNum).FilterBank_Pulses_strct(BAT(BatNum).CurrPulseNum);
     fc_detect_struct = BAT(BatNum).fc_detect_struct;
    
     integrated_pp = BAT(BatNum).integrated_detect.int_signal;
     pks_int = BAT(BatNum).integrated_detect.pks;
     locs_int = BAT(BatNum).integrated_detect.delays;
  
     tt= [1:numel(BAT(BatNum).integrated_detect.int_signal)];
 end % if main_run
 %%% 
    tt_start_time = rx_strct.TxStartTime;
    
    max_rx = max([jam(:); rx(:)]);
    k_factor = max(integrated_pp)/(10*log10(max_rx));
    
%     k_factor = 1 ; % numel(unique_active_fcs);
    
plot_db_flag= 0;
 if plot_db_flag
     sig=  10*log10(integrated_pp / k_factor);
     pks =  10*log10(pks_int/ k_factor);
 else
     sig = integrated_pp / k_factor;
     pks = pks_int/ k_factor;
 end % if plot_db
 
%  [yp,tp] = findpeaks(integrated_pp, 'MinPeakDistance',min_diff_peaks, ...
%         'MinPeakHeight', 5*NoiseLevel, 'MinPeakProminence', min_peak_ct)
 plot(tt+tt_start_time, sig,'.-b')
 plot(locs_int + tt_start_time, pks, 'dk')

 hold on
 plot(fc_detect_struct.rise_time_det, fc_detect_struct.rx_power_db , 'dr')
 xlabel('samples')
 ylabel('amp')
 
 
%  figure
%  hold on
%  plot(tt+tt_start_time,10*log10(rx_vec2),'ok')
%  plot(tt+tt_start_time,10*log10(rx_vec),'b')


 
%% response in current freq 
figure
hold on
title(num2str(k_idx))

tt_start_time = rx_strct.TxStartTime;
% total
% plot( tt+tt_start_time,10*log10(rx_vec+0.01),'.-k')
plot(tt_start_time, 0, 'dk', 'MarkerSize', 12)
plot(tt_start_time+tt(end), 0, 'dk', 'MarkerSize', 12)
plot( tt+tt_start_time,rx_vec,'-k')

%Signal only
% plot( tt+tt_start_time,10*log10(FilterBankRxMat(k_idx,tt_ind)+0.01),'.-b')
plot( tt+tt_start_time,FilterBankRxMat(k_idx,tt_ind),'.-b')
% jam only
plot( tt+tt_start_time,10*log10(JammingMat(k_idx,tt_ind)+0.01),'.-r')

% plot(amp_ax, tt, dy,'r','LineWidth',1)
% plot(amp_ax, tt(idt), dypks, 'or');
if ~isempty(rise_times)
% % %     plot( tt(idt)+tt_start_time, rx_vec(idt),'dr');
    plot( tt(rise_times)+tt_start_time, rx_vec(rise_times),'*r');
%      plot( tt(idt)+tt_start_time, rx_vec(idt),'*r');
end % if ~isempty(rise_times)

% the current fcx_times
v_fcs = 1:numel(FilterBank.filter_fc);
for k = 1:numel(rx_strct.RxEchosFromPrey)
    fcs_idx =  round(interp1(FilterBank.filter_fc, v_fcs, rx_strct.RxEchosFromPrey(k).Rx_freqs))';
    hold on
    plot( rx_strct.RxEchosFromPrey(k).Rx_times, fcs_idx/10, 'g')
end % for k

% the integrated estimation
plot( tt+tt_start_time , integrated_detect.int_signal,'g')
plot( tt(integrated_detect.delays)+tt_start_time, integrated_detect.pks, '*g') 

%% the peaks
figure ;

findpeaks(rx_vec, 'MinPeakDistance',min_diff_peaks, ...
'MinPeakHeight', NoiseLevel, 'MinPeakProminence', min_peak_ct, 'Annotate','extents')
title(num2str(k_idx))
 [dypks,idt, wdths, proms] = findpeaks(rx_vec, 'MinPeakDistance',min_diff_peaks, ...
                  'MinPeakHeight', NoiseLevel, 'MinPeakProminence', min_peak_ct)



%% detection
ind = rx_vec > NoiseLevel;
% plot(, ones(size(tt(ind))) * fc_detect_struct(kk).active_fc/1000, '*')

% plot(tt, (rx_vec > NoiseLevel) * fc_detect_struct(kk).active_fc/1000, 'r*-')
% Tx freq
plot(amp_ax, rx_strct.TxTimes_us, rx_strct.TxActiveFcs/1000)
uTx_Fcs_idx = unique(rx_strct.TxActiveFcs_idx,'stable')';
uTx_Fcs = unique(rx_strct.TxActiveFcs,'stable')'/1000;
tx_fcs_st = rx_strct.TxTimes_us(find(rx_strct.TxActiveFcs_idx == uTx_Fcs_idx(kk), 1))
plot(amp_ax, tx_fcs_st, uTx_Fcs(kk), '*') 

 

%% delay plot for decision % run in ecoes_decision

u_all_delays = unique(all_delays);
all_elays_counts = nonzeros(accumarray(all_delays',all_delays,[],@numel));

figure 
hold on
title(['delay tollerance th = ', num2str(delay_tollerance)])
plot(detected_delays, detected_counts(detections_idx), 'rd', 'MarkerSize', 8)

bar(detected_delays_temp, detected_counts,'r')
plot(u_all_delays, all_elays_counts,'bo')
plot(detected_echos_st_ref, fb_detect_th*ones(size(detected_echos_st_ref)), 'go', 'MarkerSize', 10)

% bar(u_all_delays, all_elays_counts,'b')

% plot(detected_delays_temp, detected_counts,'r*')
plot([u_all_delays(1),u_all_delays(end)], ones(1,2)*fb_detect_th,'r--')

L = legend('detected','clustered', 'all', 'power det');


% the delay detected by integrator
plot(integrated_detect.delays, fb_detect_th* ones(size(integrated_detect.delays)), 'g*', 'MarkerSize', 10)

%% Signals , Jamming and detections after FilterBank_detector and echoes_decision

% CurrentPulseNum = BAT(BatNum).CurrPulseNum;
% curr_start_time = BAT(BatNum).curr_start_time;
% curr_end_time = BAT(BatNum).curr_end_time;
% fc_jamming_struct = BAT(BatNum).fc_jamming_struct; 
% detected_echoes = BAT(BatNum).FilterBank_all_detecetd_echoes(CurrentPulseNum);

rx_strct = BAT(BatNum).FilterBank_Pulses_strct(CurrentPulseNum);
tt = curr_start_time: curr_end_time;
% detections
det_rise_time = [fc_jamming_struct.rise_time_det];
det_rx_db = [fc_jamming_struct.rx_power_db];
raw_delays = [fc_jamming_struct.raw_detections];

% ecoes dtected 
decided_times = [detected_echoes.all_detected_delays]+ curr_start_time;
decided_powers_db = detected_echoes.powers_db;
% echoes rx
all_rx_times = [rx_strct.RxEchosFromPrey.Rx_times];
all_fcs_idx = repmat([rx_strct.TxActiveFcs_idx]', 1, AllParams.SimParams.TotalPreysNumber);

figure
title(['detector input vs outputs, BAT, pulse:  ', num2str(BatNum), ', ',num2str(BAT(BatNum).CurrPulseNum)])
det_ax = gca;
hold on
plot(tt, 10*log10(abs(rel_jamming_mat+0.01)),'.r')
plot(tt, 10*log10(abs(rel_rx_mat+0.01)),'.k')
plot(det_rise_time, det_rx_db,'ob')
plot(det_ax, det_ax.XLim, zeros(size(det_ax.XLim)), 'r--') % reference TH

% freqs
% plot( all_rx_times, all_fcs_idx, 'g.')

bar(decided_times, detected_echoes.cluster_counts,'FaceColor','g', 'BarWidth', 0.5) % detected_echoes.cluster_width)
histogram(det_rise_time, numel(det_rise_time))

% final decision
plot(decided_times, decided_powers_db,'*r', 'MarkerSize',11, 'LineWidth',2)

% integratio filter
hold on
plot(tt, 10*log10(abs(integrated_detect.int_signal+0.01)), 'g.')
plot(tt(integrated_detect.delays), 10*log10(abs(integrated_detect.pks)),'*g')


%% check vs InterferenceFullStruct
InterferenceFullStruct = BAT(BatNum).InterferenceFullStruct;

notempty_idx = find(~cellfun(@isempty,{InterferenceFullStruct.BatTxPulseNum}) );

all_us_times = vertcat(InterferenceFullStruct(notempty_idx).Times_precise)*filterbank_fs_ts;
all_batTx = vertcat(InterferenceFullStruct(notempty_idx).BATTxNum);
all_jam_rx_powers = 10*log10(abs(vertcat(InterferenceFullStruct(notempty_idx).Power)+0.1));
all_freqs = vertcat(InterferenceFullStruct(notempty_idx).Freqs);
figure
plot(all_us_times, all_jam_rx_powers, 'g.')
% % tt2 =  BAT(BatNum).curr_start_time: (BAT(BatNum).curr_start_time + size(BAT(BatNum).FilterBank_jam_raw_update_mat,2)-1) ;
% % plot(tt2,10*log10(abs(BAT(BatNum).FilterBank_jam_raw_update_mat+0.1)),'.b')


 %% freq plot

freq_flag = 0;
if freq_flag
    figure
    hold on
    freq_ax =gca;
    title(' freqs')
    for k=1:AllParams.SimParams.TotalPreysNumber
        plot(CurrEchosFromPreyStruct.EchoDetailed(k).EchosnTimesFull*125,...
            CurrEchosFromPreyStruct.EchoDetailed(k).EchosFreqsFull ,'o-')
        plot(rx_strct.RxEchosFromPrey(k).Rx_times, rx_strct.RxEchosFromPrey(k).Rx_freqs/1000, '.')
    end
    plot(tt, (rx_vec > NoiseLevel) * fc_detect_struct(kk).active_fc/1000, 'r*-')
    plot(tt, (rx_vect > 0) * fc_detect_struct(kk).active_fc, 'm.')
end % if freq_flag

% figure()
% hold on
% title('Delays')
% plot(unq_delays, delays_numel,'*k')
% plot(actual_delays, ones(size(actual_delays)),'ro')

% figure
% hold on
% for nn= 1:numel(fc_detect_struct)
%     plot(fc_detect_struct(nn).detection_delays,'o')
% end %for nn= 1:numel(fc_detect_struct)