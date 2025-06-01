% active channels
active_channels = ~isnan([FilterBank_Out.Ref_ipks]);
nch = sum(active_channels);

%% The Detector
% all_detections

all_ipks = [FilterBank_Out.Channels.ipks_shift];

% histogram of detections
bin_size = 1; % 50e-6*fs;% 0.1/fs1e-3*fs;
[N,edges] = histcounts(all_ipks,'BinWidth',bin_size);
if bin_size == 1
    mid_bins = edges(1:end-1);
else 
    mid_bins = edges(1:end-1) + diff(edges)/2;
end % if
% Gaussian Filter
L =64;
sigma = 50e-6*fs; % the resolutin of detections 50 microsec
alpha = (L-1)/(2*sigma);
g_win = gausswin(L,alpha);

% Apply the filter
smoothN = conv(N, g_win,"same");

%Final Detections
pkTH =  nch/3;
[pks, ipks] = findpeaks(smoothN, 'MinPeakDistance', sigma,  'MinPeakHeight', pkTH);
detetcts_times = edges(ipks)/fs

%% Figure
figure
hold on
bar(mid_bins/fs,N)
plot(mid_bins/fs,smoothN,'r')
plot(mid_bins(ipks)/fs, pks, '*r', 'MarkerSize',8)
% THersholds
plot(xlim, nch/2*[1,1] , 'r--')
plot(xlim, nch/3*[1,1] , 'r--')

plot(xlim, nch/4*[1,1] , 'r--')


title('20dB diff 0.5mec')
xlabel('time (msec)')
ylabel('number of detections')

%%
50 mirc
%% 
bin_size =1;
tic

for k =1:1000
    [N,edges] = histcounts(all_ipks,'BinWidth',bin_size);
    NN = conv(N, w,"same");
end

toc
