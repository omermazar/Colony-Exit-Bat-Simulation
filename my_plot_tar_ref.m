function [] = my_plot_tar_ref(mat_target, FR_target, ix_start_target, f, fs)

%% target refernce
t = (1 : size(mat_target,1) ) ;
% move to start of the signal and msec
t = (t - ix_start_target ) / fs * 1000; 
figure
% spectrogram
subplot(2,1,1)
imagesc(t, round(f), mat_target')
set(gca, 'YDir', 'normal')
xlabel('time(ms)')
ylabel('freq(Hz)')
% hold on
% plot([nstart-nresponse, nstart+nresponse], [min(f), max(f)], 'w--', 'LineWidth', 10)
title('PS: Response')
% time- response
% subplot(2,1,2)
% hold on
% plot(t(nstart-1 + (1: (numel(sig_target)))), sig_target)
% plot(t, sum(mat_target')./ max(sum(mat_target'),[],'all' ),'k', 'LineWidth', 2 )
% title('Shifted Time Response')
% frequency response
subplot(2,1,2)
plot(round(f), FR_target, 'r', 'LineWidth', 1.5)
title('FR - freqeuncy Response')
sgtitle('the target reference')
xlabel('freq(Hz)')
ylabel('Response')
