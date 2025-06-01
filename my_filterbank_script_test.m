

%% Filter-Bank Params
num_of_filters = 120; %80;
fs = 200e3;
f_low = 20e3;
f_high = 80e3;
response_time = 0.5e-3; % good for f_low = 20kHz
fft_resolution = 1e3;
tir = 0:1/fs:response_time;

% LPF 
f_lowpass= 8e3; % 8e3;
[LPF_filt.b_filt, LPF_filt.a_filt] = butter(6, f_lowpass/(fs/2)); % butter(6, f_lowpass/(fs/2));



%% create the filter bank 
[AllFiltersParams]  = FilterBank_Acous_Resp(num_of_filters, fs, f_low, f_high, response_time, fft_resolution)
fc  = AllFiltersParams.filter_fc;


%% input signals
% nfilt = 10;
% fck = fc(nfilt);


call_dur = 7e-3; % 5e-3; %
PulsePower = 1;
tin = [0:1/fs:call_dur]; % time vector
nwin = 256;
fchirp_high = 47e3; % 60e3;
fchirp_low = 40e3;

chirpCall = PulsePower * chirp(tin, fchirp_high, max(tin), fchirp_low, 'logarithmic'); % chirp(tin, f_high, max(tin), f_low, 'logarithmic')
% chirpCall = PulsePower * chirp(tin, f_high, max(tin), f_low, 'linear');

% chirp_freqs = f_low*(f_high./f_low).^(tin/call_dur);
chirp_freqs = fchirp_high*(fchirp_low./ fchirp_high).^(tin/call_dur);

chirpCall_zer = [zeros(1,nwin), chirpCall ,zeros(1,nwin)]; % for the plot
% spectrogram(chirpCall_zer, nwin, nwin-1, nwin, fs,'yaxis', 'MinThreshold',-50)
[~,f,tt,ps] = spectrogram(chirpCall, nwin, nwin-1, nwin, fs,'yaxis', 'MinThreshold',-50); % for the figure


% sig_sin = sin(2*pi*fck*tin);
% sig_sin_off = sin(2*pi*1.5*fck*tin);
sig_in = chirpCall; %sig_sin_off; % sig_sin; %chirpCall;

%% Single echo
ipi = 50e-3;
nlvl = 1/sqrt(1000); %0.01;
t_ipi = 0:1/fs:ipi;
sig_out = nlvl*randn(size(t_ipi));
t_echo = 10e-3;
ix = t_echo*fs:(t_echo*fs+numel(sig_in)-1);
sig_out(ix) = sig_out(ix)+sig_in;

%% Multi Rx sig
ipi = 0.05;
t_ipi = 0:1/fs:ipi; %3 ipi
nlvl = 1/sqrt(100); %0.01;
sig_out_m = nlvl*randn(size(t_ipi));
t_diff = 2e-3;
start_time = 10e-3;

for k=1:2
    ix2 = round(start_time*fs + ((k-1)*t_diff*fs:((k-1)*t_diff*fs+numel(sig_in)-1)));
    sig_out_m(ix2) = sig_out_m(ix2) + 0.5*sig_in; % 10^(-(k-1))*sig_in; % sig_in; % 10^(-(k-1))*sig_in; % 0.5^k * sig_in; 
%     sig_out_m(ix2) = sig_out_m(ix2)+0.5^k * sig_in;
end % for k

%% Strong Signal 
% sig_out_m = [sig_out_m, nlvl*randn(size(t_ipi))];
sig_out_m = sig_out;
start_time = t_echo + 7e-3;
ix3 = round(start_time*fs : (start_time*fs + numel(sig_in)-1) );
sig_out_m(ix3) = sig_out_m(ix3) + sqrt(10^4)*sig_in;  % 30% sig_in; % 10^(-(k-1))*sig_in; % 0.5^k * sig_in; 

%% Interference in other freqs
% sig_out_m = nlvl*randn(size(sig_out));
sig_out_m = sig_out;

time_diff = 1e-3;%  2.5e-3;
start_time = t_echo + time_diff;
fint_low = 45e3; 
fint_high = 52e3;
int_lvl = 10; % dB
IntCall =  PulsePower * chirp(tin, fint_high, max(tin), fint_low, 'logarithmic'); % chirp(tin, f_high, max(tin), f_low, 'logarithmic')

ix3 = round(start_time*fs : (start_time*fs + numel(sig_in)-1) );
sig_out_m(ix3) = sig_out_m(ix3) + 10^(int_lvl/20)*IntCall;  % 30% sig_in; % 10^(-(k-1))*sig_in; % 0.5^k * sig_in; 


%% OUTPUT
[FilterBank_Ref] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_in, fs,1)
[FilterBank_Out] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_out_m, fs, false, FilterBank_Ref)

% % %% output - one channel - input
% % nfilt = 40;
% % fck = fc(nfilt);
% % 
% % out1= conv(sig_in, AllFiltersParams.filter_impulse_resp(nfilt,:) );
% % tout = ( 1:numel(out1) ) / fs; %sec
% % % rectifier
% % out2 = zeros(size(out1));
% % ix_pos = out1 > 0;
% % out2 = zeros(size(out1));
% % out2(ix_pos) = out1(ix_pos);
% % % LPF 
% % out3 = filter(b_filt, a_filt, out2);
% % 
% % % the time estimation of the freq in 
% % t_fcs = interp1(chirp_freqs, tin, fc);
% % [pks, ipks] = findpeaks(out3, 'MinPeakDistance',20,  'MinPeakHeight',0.5, 'NPeaks', 1);
% % delay = tout(ipks);
% % 
% % %% plot the channel
% % figure
% % ax1= subplot(2,1,1);
% % % spectrogram(chirpCall_zer, nwin, nwin-1, nwin, fs,'yaxis', 'MinThreshold',-50)
% % imagesc(ax1, tt*1e3, f, flipud( ps) )
% % % plot(ax1, tin*1e3, sig_in);
% % title('Input')
% % 
% % ax2 = subplot(2,1,2); hold on
% % plot(ax2, tout*1e3, out1,'DisplayName', ['BPF:' num2str(round(fck))])
% % plot(ax2, tout*1e3, out2,'DisplayName','Rectifier')
% % plot(ax2, tout*1e3, out3, 'k', 'LineWidth', 1.5, 'DisplayName','LPF')
% % xlabel('time (msec)')
% % legend
% % title(['Output, channel:' num2str(nfilt)])
% % 
% % linkaxes([ax1, ax2, ax3],'x')
% % 
% % %% plot several channels
% % figure
% % ax1= subplot(2,1,1);
% % % spectrogram(chirpCall_zer, nwin, nwin-1, nwin, fs,'yaxis', 'MinThreshold',-50)
% % imagesc(ax1, tt*1e3, f,  ps )
% % set(ax1,'YDir','normal')
% % 
% % 
% % title('Input')
% % 
% % ax2 = subplot(2,1,2); hold on
% % title(ax2, 'Output channels')
% % nfilt = (5:10:max(numel(fc)));
% % for ifilt = nfilt
% %     fck = fc(ifilt);
% %     % the response
% %     out1= conv(sig_in, AllFiltersParams.filter_impulse_resp(ifilt,:) );
% %     tout = ( 1:numel(out1) ) / fs; %sec
% %     out2 = zeros(size(out1));
% %     ix_pos = out1 > 0;
% %     out2(ix_pos) = out1(ix_pos);
% %     out3 = filter(b_filt, a_filt, out2);
% %     plot(ax2, tout, out2, 'r');
% %     plot(ax2, tout, out3, 'LineWidth', 1.5, 'DisplayName',['Ch:', num2str(ifilt), ', freq:~', num2str(round(fck/1000))]);
% % end % for igilt
% % Leg = legend;
% % %Leg.String = Leg.String(2:2:end);
% % linkaxes([ax1, ax2],'x')


% % %% All channels - Sig in
% % figure
% % ax1= subplot(2,1,1);
% % % spectrogram(chirpCall_zer, nwin, nwin-1, nwin, fs,'yaxis', 'MinThreshold',-50)
% % [~,f,tt,ps] = spectrogram(sig_in, nwin, nwin-1, nwin, fs,'yaxis', 'MinThreshold',-50); % for the figure
% % 
% % imagesc(ax1, tt, f,  ps )
% % set(ax1,'YDir','normal')
% % ax1.XLabel.String = 'time(msec)';
% % ax1.YLabel.String = 'Freq (kHz)';
% % ax1.YTick = [0: 10e3: max(fc)+10e3];
% % ax1.YTickLabel = ax1.YTick/1000;
% % grid on
% % title(ax1, 'input')
% % 
% % ax2 = subplot(2,1,2); hold on
% % nchannels = numel(fc);
% % resp_time = nan(1,nchannels);
% % resp_pks = nan(1,nchannels);
% % for ifilt = 1:nchannels
% %     fck = fc(ifilt);
% %     % the response
% %     out1= conv(sig_in, AllFiltersParams.filter_impulse_resp(ifilt,:) );
% %     tout = ( 1:numel(out1) ) / fs; %sec
% % 
% %     out2 = zeros(size(out1));
% %     ix_pos = out1 > 0;
% %     out2(ix_pos) = out1(ix_pos);
% %     out3 = filter(b_filt, a_filt, out2);
% %     minPeakDist = 2.5e-4*fs;
% %     [pks, ipks] = findpeaks(out3, 'MinPeakDistance',minPeakDist,  'MinPeakHeight',0.5); % , 'NPeaks', 1);
% %     if ~isempty(pks) % else is nan
% %         [mpk, impk] = max(pks); % the highest peak
% %         resp_time(ifilt) = ipks(impk); % 1/fs*
% %         resp_pks(ifilt) = pks(impk);
% %     end % if ~isempty(pks) 
% % %     plot(ax2, tout, out2, 'r');
% %     plot(ax2, tout, out3, 'LineWidth', 1.5, 'DisplayName',['Ch:', num2str(ifilt), ', freq:~', num2str(round(fck/1000))]);
% %     plot(ax2, 1/fs*resp_time(ifilt), resp_pks(ifilt), 'or');
% % end % for igilt
% % linkaxes([ax1, ax2],'x')
% % 
% % %%% plot the peaks as function pf channel
% % plot(fc,resp_pks,'*k-', 'DisplayName', [num2str(call_dur*1e3),'ms'])
% % title('Chirp 40-20kHz')
% % legend
% % %% sum of all channels
% % figure
% % ax1= subplot(2,1,1);
% % % spectrogram(chirpCall_zer, nwin, nwin-1, nwin, fs,'yaxis', 'MinThreshold',-50)
% % [~,f,tt,ps] = spectrogram(sig_out, nwin, nwin-1, nwin, fs,'yaxis', 'MinThreshold',-100); % for the figure
% % 
% % imagesc(ax1, tt, f,  ps )
% % set(ax1,'YDir','normal')
% % ax1.XLabel.String = 'time(msec)';
% % ax1.YLabel.String = 'Freq (kHz)';
% % ax1.YTick = [0: 10e3: max(fc)+10e3];
% % ax1.YTickLabel = ax1.YTick/1000;
% % grid on
% % title(ax1, 'input')
% % 
% % ax2 = subplot(2,1,2); hold on
% % nchannels = numel(fc);
% % % % % resp_time = zeros(1,nchannels);
% % % % % resp_pks = zeros(1,nchannels);
% % out(nchannels) = struct('response', [], ...
% %     'shifted_response',[], ...
% %     'pks',[], ...
% %     'ipks',[]);
% % 
% % 
% % for ifilt = 1:nchannels
% %     fck = fc(ifilt);
% %     % the response
% %     out1= conv(sig_out, AllFiltersParams.filter_impulse_resp(ifilt,:) );
% %     tout = ( 1:numel(out1) ) / fs; %sec
% % 
% %     out2 = zeros(size(out1));
% %     ix_pos = out1 > 0;
% %     out2(ix_pos) = out1(ix_pos);
% %     out3 = filter(b_filt, a_filt, out2);
% %     minPeakDist = 2.5e-4*fs;
% %     [pks, ipks] = findpeaks(out3, 'MinPeakDistance',minPeakDist,  'MinPeakHeight',0.5); % , 'NPeaks', 1);
% % %     [mpk, impk] = max(pks); % the highest peak
% % %     resp_time(ifilt) = 1/fs*ipks(impk);
% % %     resp_time(ifilt) = ipks;    
% % %     resp_pks(ifilt) = pks; % (impk);
% % %     plot(ax2, tout, out2, 'r');
% %     plot(ax2, tout, out3, 'LineWidth', 1.5, 'DisplayName',['Ch:', num2str(ifilt), ', freq:~', num2str(round(fck/1000))]);
% %     plot(ax2, 1/fs*ipks, pks, 'or');
% %     
% %     % output
% %     out(ifilt).response = out3;
% %     out(ifilt).shifted_response = [out3(resp_time(ifilt) : end),  zeros(1,resp_time(ifilt)-1) ];
% %     out(ifilt).pks = pks;
% %     out(ifilt).ipks = ipks;
% %     plot(ax2, tout, out(ifilt).shifted_response, '--', 'LineWidth', 1.5, 'DisplayName',['Ch:', num2str(ifilt), ', freq:~', num2str(round(fck/1000))]);
% % end % for igilt
% % linkaxes([ax1, ax2],'x')
% % 
% % % sum
% % sr = sum(vertcat(out.shifted_response));
% % sr = sr/max(sr); % normalize
% % figure
% % plot(tout, sr, '.-') 
% % hold on
% % [pks_sum, ipks_sum, ] = findpeaks(sr, 'MinPeakDistance',minPeakDist,  'MinPeakHeight',0.5); % , 'NPeaks', 1);
% % plot(ipks_sum/fs, pks_sum, '*')
% % %% test time
% % 
% % tic
% % for k=1:100
% %     [FilterBank_Out_m] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_out_m, fs, false, FilterBank_Ref);
% % end
% % toc