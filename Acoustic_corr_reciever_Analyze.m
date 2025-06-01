%% Test Correlation Receiver
fig = figure;

warning('off', 'signal:findpeaks:largeMinPeakHeight')
warning('off', 'MATLAB:catenate:DimensionMismatch')
warning('off', 'MATLAB:inpolygon:ModelingWorldLower')
warning('off', 'MATLAB:polyshape:repairedBySimplify')
warning('off', 'MATLAB:polyshape:boolOperationFailed')
warning('off', 'signal:findpeaks:largeMinPeakHeight')
warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode')
warning('off', 'MATLAB:legend:IgnoringExtraEntries')

ax1 = subplot(3,1,1) ; hold on
ax2 = subplot(3,1,2) ; hold on; grid minor
ax3 = subplot(3,1,3); hold on;

fig.Position = [ 50 150 1500 600];
ax1.Position =  [0.13   0.55    0.7750    0.38];
set(ax2,'Position', [0.13   0.35    0.7750    0.15]); % ax2.Position =  [0.13   0.35    0.7750    0.15];
ax3.Position =  [0.13   0.11    0.7750    0.18];

Ts = 1e3/FsAcoustic; %%% time(msec)
FilterBank_FlagXX = strcmp(AllParams.BatSonarParams.ReceiverType, 'FilterBank');

if FilterBank_FlagXX
    txt_sup = 'FB';
    timesW =  0:(numel(FilterBank_RxWanted.Sum_out)-1);
    timesA =  0:(numel(FilterBank_RxAll.Sum_out)-1);
%     int_times = (0:(numel(int_cors)-1))*Ts; % times_int; % FilterBank_RxInt.time_vec*1000;
else
    txt_sup = 'Corr';
    timesW = 0:(numel(corsW)-1);
    timesA = 0:(numel(corsA)-1);
%     int_times = 0:(numel(int_cors)-1);
end % if FilterBank_FlagXX
%% Set colors
[colStruct] = myColorScheme(AllParams.SimParams.TestMode);

%%% Orignial Masking
orig_nTimes = BAT(BatNumRx).FindsMaskingStruct(PulseNum).CoverTimeJamTimes;
orig_MaskedPreys =  BAT(BatNumRx).FindsMaskingStruct(PulseNum).MaskedPreys;

sgtitle([txt_sup, ', Bat:', num2str(BatNumRx), '; PulseNum:', num2str(PulseNum), ...
    '; ', BAT(BatNumRx).ManueverCmdStruct(PulseNum).ManueverStage, ...
    ',  Orig:', num2str(numel(DetectedEchoesVec)), '\' , num2str(numel(orig_MaskedPreys)), ...
    ',  Acous:', num2str(numel(wanted_detected_times)), '\' , num2str(sum(jamm_flg))]) 

%% All &  INTERFERENCE

pre_detected_rxlvl = 10*log10(abs(CurrEchosStruct.EchosAttenuation(ix_detections)))+PulsePower;
norm_fact1 = max([10*log10(CurrEchosStruct.EchosAttenuation)+PulsePower, pre_detected_rxlvl]) ;
if ~FilterBank_FlagXX
    switch detectMode
        case 'envelope'
            corrAll    = envelope(corsA, 64, 'analytic');
%             intCors    = envelope(int_cors, 64, 'analytic');
            wantedCors = envelope(corsW, 64, 'analytic');
        case 'none'
            corrAll    = corsA;
%             intCors    = int_cors;
            wantedCors = corsW;
    end % switch

    plot(ax1, timesA*Ts, corrAll, 'k', 'LineWidth' ,1, 'DisplayName', 'Correlation All')
%     plot(ax1, int_times*Ts, intCors, 'r--', 'LineWidth' ,1.5 , 'DisplayName', 'Correlation Interf')
else % if ~FilterBank_FlagXX
    plot(ax1, FilterBank_RxAll.time_vec*1000, FilterBank_RxAll.Sum_ipks_hist, 'k', 'LineWidth' ,1, 'DisplayName', 'Correlation All') % /norm_fact1
%     plot(ax1, FilterBank_RxInt.time_vec*1000, FilterBank_RxInt.Sum_ipks_hist, 'r--', 'LineWidth' ,1.5 , 'DisplayName', 'Correlation Interf')
end % if ~FilterBank_FlagXX


%% Wanted
% plot( (1:numel(Rx_Acoustic_WantedSig))*Ts, Rx_Acoustic_WantedSig,'b','DisplayName', 'Wanted')
% plot( (1:numel(Tx_Acoustic_Call))*Ts, Tx_Acoustic_Call,'g', 'DisplayName', 'Tx')

try 
    plot(ax1, pre_detect_times*Ts, pre_detected_rxlvl*0.001,'*m', 'DisplayName', ' all pre-detections')
catch
    disp('Error: did not plot pre_detected_rxlvl')
end
try
    plot(ax1, pre_detect_times(~detected_ix)*Ts, pre_detected_rxlvl(~detected_ix)*0.001,'sm',...
        'MarkerSize',10, 'LineWidth', 3, 'DisplayName', 'Missed by Corr')
catch
    disp('Error: did not plot Missed by Corr')
end
try
    plot(ax1, CurrEchosStruct.EchosTimes*FsAcoustic*SampleTime*Ts, (10*log10(CurrEchosStruct.EchosAttenuation)+PulsePower)*0.001, 'om', 'DisplayName', 'All Echoes')
    catch
    disp('Error: did not plot All Echoes')
end


% % plot(wanted_pks_times*Ts, wanted_pks./norm_fact1, 'go', 'DisplayName', 'Detected wanted')
% % plot(wanted_detected_times*Ts, wanted_detected_pks./norm_fact1, 'og', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'matching')

if ~FilterBank_FlagXX
    plot(ax1, wanted_pks_times*Ts, wanted_pks, 'go', 'DisplayName', 'Detected wanted')
    plot(ax1, wanted_detected_times*Ts, wanted_detected_pks, 'og', 'MarkerFaceColor', 'g', 'MarkerSize', 7, 'DisplayName', 'matching')
    plot(ax1, timesW*Ts, wantedCors, 'g--', 'LineWidth' ,1.5, 'DisplayName', 'Correlation Wanted')
    
    
else % if ~FilterBank_FlagXX
    plot(ax1, wanted_pks_times*Ts, wanted_pks, 'go', 'DisplayName', 'Detected wanted')
    plot(ax1, wanted_detected_times*Ts, wanted_detected_pks, 'ok', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'matching')
    plot(ax1, all_detected_times/FsAcoustic*1000, all_detected_pks, 'ks', 'MarkerFaceColor', 'g', 'DisplayName', 'Estimated detect')
    plot(ax1, int_Times/FsAcoustic*1000, int_Lvl, 'dr', 'MarkerFaceColor', 'r', 'DisplayName', 'Masking peak')
    plot(ax1, FilterBank_RxWanted.time_vec*1000, FilterBank_RxWanted.Sum_ipks_hist, 'g--', 'LineWidth' ,1, 'DisplayName', 'Correlation Wanted')
    active_channels = ~isnan([FilterBank_RxWanted.Ref_ipks]);
    nch = sum(active_channels);
    corrTH = FilterBank_RxWanted.detection_TH; % /3;
end % if ~FilterBank_FlagXX

% Rx Detection TH
plot(ax1, [min(timesA), max(timesA)]*Ts, ones(1,2)*DetectTH, 'm:', 'LineWidth', 1.5, 'DisplayName', 'RxLvl TH')
% Corr detection TH
plot(ax1, [min(timesA), max(timesA)]*Ts, ones(1,2)*corrTH, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Corr TH')
ax1.YLim(1) = -1;
% Legend
leg1 = legend(ax1);
% xlabel('time (msec)')
ylabel(ax1, 'NormDetector  / RxdB')

% text
xLim = xlim; yLim = ylim;
x1 = ax1.XLim(1) + 0.05*diff(ax1.XLim);
y1 = ax1.YLim(2) - 0.05*diff(ax1.YLim);
y2 = y1 - 0.06*diff(ax1.YLim);
text(ax1, x1, y1, [num2str(pre_detect_times*Ts, 3), 'ms'], 'FontSize', 7);
text(ax1, x1, y2, [num2str(pre_detected_rxlvl, 3), 'dB'], 'FontSize', 7);


%% Masking
% plot(all_pks_times*Ts, all_pks, 'ks', 'MarkerSize',8, 'LineWidth',2 , 'DisplayName', 'Detected all')
plot(ax1, direct_jammed_times*Ts, direct_jammed_pks, '+r', 'MarkerSize', 12, 'LineWidth',3, 'DisplayName', 'Direct Jam')
plot(ax1, classifier_jammed_times*Ts, classifier_jammed_pks, '*r', 'MarkerSize', 6, 'LineWidth',2, 'DisplayName', 'Classifier Jam')
plot(ax1, wanted_detected_times(jamm_flg)*Ts, wanted_detected_pks(jamm_flg)./norm_fact1 , 'or', 'MarkerSize', 12, 'LineWidth',3, 'DisplayName', 'FwdBwd Jam')
%plot(masking_pks_times*Ts, masking_pks_lvls, 'ro', 'MarkerSize',4, 'LineWidth',1 , 'DisplayName', 'all masking')

%% CAlssifier - Detections and False Alarms

plot(ax1, DetectedCalissified_nTimes*Ts, zeros(size(DetectedCalissified_nTimes)), '+b', 'MarkerSize', 9, 'LineWidth',2, 'DisplayName', 'ClassifierDetection')

plot(ax1, FalseAlarmsClassified_nTimes*Ts, zeros(size(FalseAlarmsClassified_nTimes)), 'db', 'MarkerSize', 9, 'LineWidth',3, 'DisplayName', 'ClassifierFA')

%%% Orignial Masking
% orig_nTimes = BAT(BatNumRx).FindsMaskingStruct(PulseNum).AcousticMaskingTimes;
orig_mask_times = BAT(BatNumRx).FindsMaskingStruct(PulseNum).AcousticMaskingTimes/ FsAcoustic*1000; % (orig_nTimes-TxnTime)*SampleTime*1000; % msec
orig_mask_times = orig_mask_times -TxnTime*SampleTime*1000;

orig_mask_powers = BAT(BatNumRx).FindsMaskingStruct(PulseNum).DetectedPreysRxPower;
mask_ix = ismember( BAT(BatNumRx).FindsMaskingStruct(PulseNum).DetectedPreys, BAT(BatNumRx).FindsMaskingStruct(PulseNum).MaskedPreys);
if ~isempty(orig_mask_times)
%     orig_mask_times = orig_mask_times(mask_ix);
%     orig_mask_powers = orig_mask_powers(mask_ix);
    if FilterBank_FlagXX
        ixPowers = find(ismembertol(FilterBank_RxAll.time_vec*1000 , orig_mask_times, 0.001,'DataScale',1));
        plot (ax1, orig_mask_times,  FilterBank_RxAll.Sum_ipks_hist(ixPowers), '+m', 'MarkerSize', 10, 'LineWidth',2, 'DisplayName', 'Original Masking')
    else
        ixPowers =  find(ismembertol(timesA*Ts , orig_mask_times, 0.001,'DataScale',1));
        plot (ax1, orig_mask_times, corrAll(ixPowers), '+m', 'MarkerSize', 10, 'LineWidth',2, 'DisplayName', 'Original Masking')
    end
end % if ~isempty(orig_mask_times)

% plot(timesw*Ts, envelope(corsw,  8,'analytic'), 'b', 'LineWidth' ,1, 'DisplayName', 'Correlation')
% plot(timesw*Ts, corsw, 'r', 'LineWidth' ,1, 'DisplayName', 'Correlation')

% Interference
% plot(int_times*Ts, envelope(int_cors,64,'analytic'), 'r.', 'LineWidth' ,1, 'DisplayName', 'Correlation All')

%% The Signals
% if AllParams.SimParams.MaskingByClutter
%    norm_fact2 = max([20*log10(abs(Rx_Acoustic_Interference)) , 20*log10(abs(Rx_Acoustic_WantedSig)), 20*log10(abs(Rx_Acoustic_Clutter))]);
% else
norm_fact2 = max([20*log10(abs(Rx_Acoustic_Interference)) , 20*log10(abs(Rx_Acoustic_WantedSig))]);
% end

tt = 1:numel(Rx_Acoustic_Interference	);
plot(ax2, tt*Ts, smooth(20*log10(abs(Rx_Acoustic_Interference)+1)), 'color', colStruct.masking, 'DisplayName', 'Rx_Acoustic_Interference')
plot(ax2, tt*Ts, smooth(20*log10(abs(Rx_Acoustic_WantedSig)+1)),'color', colStruct.detections , 'LineWidth', 1,'DisplayName', 'Wanted')
% plot(ax2, (1:numel(Tx_Acoustic_Call))*Ts, 20*log10(abs(Tx_Acoustic_Call)),'g', 'DisplayName', 'TxCall')
% if AllParams.SimParams.MaskingByClutter
%     plot(ax2, tt*Ts, 20*log10(abs(Rx_Acoustic_Clutter)),'m', 'DisplayName', 'Clutter')
% end

% plot( (1:numel(Rx_Acoustic_Sig))*Ts, Rx_Acoustic_Sig,'m', 'DisplayName', 'RxAll')
% if ~FilterBank_FlagXX
%     plot(ax2, timesA*Ts, envelope(corsA), 'k--', 'LineWidth' ,1, 'DisplayName', 'Correlation All')
%     plot(ax2, int_times*Ts, envelope(int_cors), 'k.', 'LineWidth' ,1.5, 'DisplayName', 'Correlation Interference')
%     plot(ax2, timesW*Ts, envelope(corsW,  64,'analytic'), 'g.', 'LineWidth' ,1.5, 'DisplayName', 'Correlation Wanted')
% else
% %     plot(ax2, FilterBank_RxWanted.time_vec*1000, FilterBank_RxAll.Sum_ipks_hist./norm_fact2, 'k--', 'LineWidth' ,1, 'DisplayName', 'Correlation All')
%     plot(ax2, timesA*Ts, FilterBank_RxInt.Sum_out./norm_fact2, 'k.', 'LineWidth' ,1.5, 'DisplayName', 'Correlation Interference')
%     plot(ax2, timesA*Ts, FilterBank_RxWanted.Sum_out./norm_fact2, 'g.', 'LineWidth' ,1.5, 'DisplayName', 'Correlation Wanted')
% end % if ~FilterBank_FlagXX

leg2 = legend(ax2);
xlabel(ax2, 'Time (msec)')
ax2.YLim(1) = -5; 
%% spectrogram
nwin = 512;
spectrogram(Rx_Acoustic_Sig, nwin, nwin-1, nwin, FsAcoustic,'yaxis', 'MinThreshold',-50);
ylim(ax3, [-0.2 80])
%% 
linkaxes([ax1, ax2, ax3], 'x')
% fig.Position = [ 50 150 1500 600];
% ax1.Position =  [0.13   0.55    0.7750    0.38];
% set(ax2,'Position', [0.13   0.35    0.7750    0.15]); % ax2.Position =  [0.13   0.35    0.7750    0.15];
% ax3.Position =  [0.13   0.11    0.7750    0.18];

