function [] = myFilterBankResponsePlot( RxSig, FilterBank_Rx, FilterBank, fs, BatDATA, BatNum, PulseNum)

activeCahnnels = find(FilterBank_Rx.Active_Channels);
freqs = FilterBank.filter_fc(activeCahnnels)/1e3;

fig = figure;
set(fig,'Units','centimeters')
fig.Position = [12   2.5  9.3  15];

%% RxLevel time plot
ax1       = subplot(4,1,1);
BAT       = BatDATA.BAT(BatNum);
startTime = BAT.TransmittedPulsesStruct(PulseNum).StartPulseTime/ 2000;
endTime   = startTime + BAT.TransmittedPulsesStruct(PulseNum).IPItoNextPulse / 2000;
ax1.YLabel.String = 'Level (dB)';
ax1.FontSize = 7;
ax1.XTickLabel = '';
My_AcousticSig_VecTest (BatDATA, BatNum, ax1, 'Obs', startTime, endTime, false);
title('')

%% Spectrogram
ax2 = subplot(4,1,2);
nwin = 512;
[~,f,t,pL] = spectrogram(RxSig, nwin, nwin-1, nwin, fs,'yaxis', 'MinThreshold',-50);
contourf(t, f/1000, 20*log10(pL+1e-5),'LevelStep', 5 , 'LineStyle','none');
ylabel('Freq (kHz)')
ax2.FontSize = 7;
ax2.YLim = [min(freqs)-3, max(freqs)+3];
ax2.XTickLabel = '';
grid on

%% Channels
ax3 = subplot(4,1,3);
hold on; grid on
for k = 1:numel(activeCahnnels)
    curChannel = activeCahnnels(k);
    detIx = find(FilterBank_Rx.Channels(curChannel).pks_shift > FilterBank_Rx.Channels(curChannel).detectionTH);
    plot(FilterBank_Rx.Channels(curChannel).ipks_shift(detIx)/fs, freqs(k)*ones(size(detIx)), '.k', 'MarkerSize',6)   %, 'MarkerFaceColor', 'k')

%     scatter(FilterBank_Rx.Channels(curChannel).ipks_shift(detIx)/fs, ...
%         freqs(k)*ones(size(detIx)), [], ...
%         FilterBank_Rx.Channels(curChannel).shifted_response(detIx), 'filled', 's') %, 'MarkerFaceColor', 'k')
end
ylabel(ax3, 'Channel Fc (kHz)')
ax3.FontSize = 7;
ax3.YLim = [min(freqs)-3, max(freqs)+3];  
ax3.XTickLabel = '';

%% FilterBank Response
ax4 = subplot(4,1,4);
grid on; hold on
plot(FilterBank_Rx.time_vec, FilterBank_Rx.Sum_ipks_hist, 'k')
plot(xlim, FilterBank_Rx.detection_TH*ones(1,2), 'r--')

xlabel(ax4, 'Time (sec)')
ylabel('Summary Response')
ax4.FontSize = 7;

%% Axis size
minY = 0.055; maxY = 0.95; 
ySize = 0.21; gap = 0.03;

linkaxes([ax1, ax2, ax3, ax4],'x')
allAx = [ax4, ax3, ax2, ax1];
for k = 1:numel(allAx)
    currAx = allAx(k);
    currAx.Position(2) = minY + (k-1) * (ySize + gap);
    currAx.Position(4) = ySize;
end