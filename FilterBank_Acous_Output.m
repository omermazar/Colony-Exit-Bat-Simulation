function [FilterBank_Out] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_in, fs, reference_flag, FilterBank_Ref, ...
    integration_method, debug_flag, AllParams, DetectedObj)

% [FilterBank_Out] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_in, fs)
% returns a strunt of the output of a FilterBank Detector fwith the input
% signal sig.
% 
% [FilterBank_Out] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_in, fs, reference_flag)
% if reference_flag is TRUE (the default is False), the Output is the
% maximal Detected peak of the peak in each channel  of the FilterBank,
% without shifting the ouput chnnels in time. if reference_flag is false the
% shiftetd channel, and FilterBank_Ref is given the output witll also shift
% the cahhnels according to the reference input signal (chirp's slope).
%       integration_method: "summation" | "channel_detection" - 'summation' for
%               summming up all canneks and than detect peoaks; "channel_detection" for seperate detection 
%               in each cannel  , according to Sallant method
%       debug_flag: true | false - for plots during running
%       AllParams: Bulti in Sturct for all parameters in the model: used
%       only for AllParams.BatSonarParams.NoiseLeveldB (0dB default)
%       DetectedObj: "Prey" (default) | "Obs" | "Conspecific" - What is the detected object?  required for
%               clutter detection (parameters for detections)
%       
% Output:
% The reposnse is: convulution with the Impulse response, rectiny, and LPF in each chanel
%   FilterBank_Out a struct with the following sields: 
%       Channels (nFiflters struct): 'response'- the respines of each cahnnel to the signal in the chanel, 
%               'shifted_response': the response of each chnnle corrected to the refence signal shape (if referrnce_flag is FALSE)           
%               'pks'- the lvl of all detected peaks in the channel 
%               'ipks' - the index pf the detected peaks (detected_time = ipeaks*fs)
%       Ref_ipks, Ref_pks - the detected peaks of the reference signals (indices and lvl, respectively)
%       Sum_out: The Filter Bank Respose: The sum of all the shifted responses 
%       Sum_pks:  the lvl of all detected peaks in the channel (normalized)
%       Sum_ipks: The indices of the peaks 

% Sep2022 - add DetectedObj - the type of object the FilteBank is detecting

%% Input
if nargin < 10
    DetectedObj = 'Prey';
end

switch DetectedObj
    case 'Prey'
      clutterMode = false;  
    case 'Obs'
      clutterMode = true;
    case 'Conspecific'
      clutterMode = false;
end % switch

if nargin < 9
    NoiseLevel = 1; %0 dBSPL
else
    NoiseLevel = 10.^(AllParams.BatSonarParams.NoiseLeveldB/10);

end

if nargin < 8
    debug_flag = 0;
end

if nargin < 7
    integration_method = "channel_detection"; %  { "summation", "channel_detection"}
end

if nargin < 6 
    FilterBank_Ref = [];
end % if nargin

if nargin < 4
 reference_flag = false;
end % if nargin


%% Genreral
fc  = AllFiltersParams.filter_fc; % the center freqs of the filterbank

%%% Detection Parameters
ref_pk_th = 1.5; % tested on several chirp - activity in the acive chanels
if ~clutterMode
    ch_pkDist  = 2.5e-4*fs; %2.5e-4*fs; %1e-3*fs; % 2.5e-4*fs; % ~the Convlution Width
    sum_pkDist = 1.5e-4*fs; % the detection threshold for the summary channel
else
    ch_pkDist  = 0.5e-4*fs; %2.5e-4*fs; %1e-3*fs; % 2.5e-4*fs; % ~the Convlution Width
    sum_pkDist = 1.5e-4*fs; % the detection threshold for the summary channel
end % if ~clutterMode

% % sig_pk_th = 1.5; % Peak detection TH
% % sig_pk_pr = 0.7; % 0.2*sig_pk_th; % prominance
% % sum_pk_th = 0.1; % the setection threshols for the sumary channel
% % sum_pk_pr = 0.025;

%%%% init
nchannels = numel(fc);
resp_time = nan(1,nchannels);
resp_pks = nan(1,nchannels);
prc = 40; % for noise-lvl estimation

% Init output
Channels(nchannels) = struct('Ref_SLL', [], ...
    'response', [], ...
    'shifted_response', [], ...
    'pks_shift', [], ...
    'ipks_shift', []);

FilterBank_Out.Active_Channels = [];
FilterBank_Out.Channels = Channels;
FilterBank_Out.Ref_ipks = nan(1,nchannels);
FilterBank_Out.Ref_pks = nan(1,nchannels);
FilterBank_Out.Ref_SLL = nan(1,nchannels); % the SLL is deifined here as the prominance of the pk/ the max prominance

FilterBank_Out.Sum_out = [];
FilterBank_Out.det_pks = [];
FilterBank_Out.det_ipks = [];
FilterBank_Out.Sum_ipks_hist = [];

if ~isempty(FilterBank_Ref)
    ref_ipks = FilterBank_Ref.Ref_ipks;
    FilterBank_Out.Ref_ipks = FilterBank_Ref.Ref_ipks;
    FilterBank_Out.Ref_pks = FilterBank_Ref.Ref_pks;
else % if ~isempty(FilterBank_Ref)
    ref_ipks = nan(nchannels);
end % if ~isempty(FilterBank_Ref)

%% Debugging and Testing Plot
% debug_flag = 1;
if debug_flag
    nwin = 512; % spectrogram
    fig = figure;
    ax1= subplot(3,1,1);
    
    % spectrogram(chirpCall_zer, nwin, nwin-1, nwin, fs,'yaxis', 'MinThreshold',-50)
    P = nextpow2(numel(sig_in))-2;
    nwin = min(nwin, 2^P);
    [~,f,tt,ps] = spectrogram(sig_in, nwin, nwin-1, nwin, fs,'yaxis', 'MinThreshold', -70); % for the figure
    
    imagesc(ax1, tt, f,  10*log10(ps) )
    set(ax1,'YDir','normal')
    ylim(ax1, [-0.2 80e3]);
    ax1.XLabel.String = 'time(msec)';
    ax1.YLabel.String = 'Freq (kHz)';
    ax1.YTick = [0: 10e3: max(fc)+10e3];
    ax1.YTickLabel = ax1.YTick/1000;
    grid on
    title(ax1, 'input')
    
    ax2 = subplot(3,1,2); hold on
    ax3 = subplot(3,1,3);hold on
   
    fig.Position = [400 200 800 550];
    linkaxes([ax1, ax2, ax3],'x')
    ax1.XLim = [0, max(tt)+nwin/fs];
end % if debug_flag

%% Minimal Detetction Lvl 
%%%%% New 02-May-2022 JAm Moth
DR = 70; % Dynamic Range (dB)
minTHLvl = max(abs(sig_in))/ 10.^(DR/20);

%% The output of each channel
for ifilt = 1:nchannels
    fck = fc(ifilt);
    
    %%%% the response
    % Covolution
    ch_out= conv(sig_in, AllFiltersParams.filter_impulse_resp(ifilt,:) );
    % Rectifier 
    ix_neg = ch_out < 0;
    ch_out(ix_neg) = 0;
    % LPF 
    ch_out = filter(LPF_filt.b_filt, LPF_filt.a_filt, ch_out);
    
    %%% detections of each channel
    if reference_flag
        
        %%%% Call signal %%%%
%         [pks, ipks] = findpeaks(ch_out, 'MinPeakDistance',ch_pkDist, 'MinPeakHeight', ref_pk_th); % , 'NPeaks', 1);
%         debug
%         [pks, ipks, ~, p ] = findpeaks(ch_out, 'MinPeakDistance',ch_pkDist, 'MinPeakHeight', ref_pk_th); % , 'NPeaks', 1);
        [pks, ipks] = myFindPeaks(ch_out, ref_pk_th, ch_pkDist ); 
        if ~isempty(pks) % else is nan
            [mpk, impk] = max(pks); % the highest peak
            ipks = ipks(impk); % 1/fs*
            pks = mpk;
            % output
            FilterBank_Out.Channels(ifilt).Ref_SLL = nan; % p / p(impk); % Changed Sep2022
            FilterBank_Out.Ref_ipks(ifilt) = ipks;
            FilterBank_Out.Ref_pks(ifilt) = pks;
        end % if ~isempty(pks)
    
    else % if reference_flag
        
        %%%% Rx Signals %%%%
        pks = nan; ipks = nan;
        switch integration_method            
            case  "summation"
                % If working with summmation no need to caculate the peaks
                             
            case "channel_detection"
                
                if ~isnan(ref_ipks(ifilt)) % only if the channle was detecetd in the reference 
                    % calcuate the peaks in each channel
                    %%%% Detection Above Noise -Level
%                     [nlvl, nstd] = noise_lvl_est(ch_out, prc); % was replaced by the real Noise Lvel for better time-performance
                    [nlvl, nstd] = deal(NoiseLevel);
                    
                    ch_pkTH = max(5 * nlvl, minTHLvl); % max(nlvl + 3*nstd, minTHLvl); % 'MinPeakHeight'
                    ch_pkPR = max(8 * nstd,  minTHLvl); %max(2.5 * nstd,  minTHLvl); % 'MinPeakProminence'
%                     [pks, ipks ] = findpeaks(ch_out, 'MinPeakDistance',ch_pkDist,  'MinPeakHeight', ch_pkTH, 'MinPeakProminence', ch_pkPR);
                    % XXX   improve sidepeaks removal
%                     [pks, ipks,~, prs ] = findpeaks(ch_out, 'MinPeakDistance',ch_pkDist,  'MinPeakHeight', ch_pkTH, 'MinPeakProminence', ch_pkPR);
                    [pks, ipks ] = myFindPeaks(ch_out, ch_pkTH, ch_pkDist ); 
%                      [pksXX, locsXX] = myFindPeaks(ch_out,ch_pkTH, ch_pkDist ); 
%                     ix_pr= pks./prs <= 5;
%                     pks = pks(ix_pr);
%                     ipks = ipks(ix_pr);
%                     % XXX   
                end
        end % switch
        
        % Shift the channel according to the reference sig, Start at the reference peak time time and zero padding the end
        if ~isnan(ref_ipks(ifilt))
            FilterBank_Out.Channels(ifilt).shifted_response = [ch_out(ref_ipks(ifilt) : end),  zeros(1, ref_ipks(ifilt)-1) ];
            FilterBank_Out.Channels(ifilt).pks_shift = pks;
            FilterBank_Out.Channels(ifilt).ipks_shift = ipks- ref_ipks(ifilt);
            FilterBank_Out.Channels(ifilt).detectionTH = ch_pkTH;
        end % if ref_ipks
                
        
    end % if reference_flag
    
    %%%% plot each channel Response %%%%%%5
    if debug_flag &&  strcmp(integration_method, "channel_detection")
        tout = ( 1:numel(ch_out) ) / fs; %sec
%         hold(ax2, 'off');
       
        if ~isnan(ref_ipks(ifilt))
%             plot(ax2, tout, ch_out, 'LineWidth', 1, 'DisplayName',['Ch:', num2str(ifilt), ', freq:~', num2str(round(fck/1000))]);
            legend(ax2, ['Ch:', num2str(ifilt), ', freq:~', num2str(round(fck/1000))])
%             plot(ax2, tout, ch_out, 'LineWidth', 1.5, 'DisplayName',['Ch:', num2str(ifilt), ', freq:~', num2str(round(fck/1000))]);
            plot(ax2, [0,tout(end)], ch_pkTH*[1,1], 'r--')
            hold(ax2, 'on');
            plot(ax2, tout,  FilterBank_Out.Channels(ifilt).shifted_response, 'LineWidth', 1);
%             plot(ax2, ipks/fs, pks, 'or')
            plot(ax2, (ipks-ref_ipks(ifilt))/fs, pks, '*r')

%             pause
            hold(ax2, 'off');
        end % if ~isnan
        

    end % if debug_flag
    
    %%%%% OutPut
    % output
    FilterBank_Out.Channels(ifilt).response = ch_out;
%     FilterBank_Out.Channels(ifilt).pks = pks;
%     FilterBank_Out.Channels(ifilt).ipks = ipks;
end % for ifilt


%% Integration of all channels

% the summation of all channels:
tout = ( 1:numel(ch_out) ) / fs; %sec
active_channels = ~isnan([FilterBank_Out.Ref_ipks]);
nch = sum(active_channels);
all_shifted_responses = reshape([FilterBank_Out.Channels.shifted_response], [], nch);
sr = sum(all_shifted_responses,2);

if ~reference_flag
    % sr = sum(vertcat(FilterBank_Out.Channels.shifted_response)); % TOO SLOW
    switch integration_method
        case  "summation"
            % If working with summmation no need to caculate the peaks
%             sr = zeros(size(ch_out));
%             nix = find(~isnan(ref_ipks));
%             for ii = nix
%                 sr = sr + FilterBank_Out.Channels(ii).shifted_response;
%             end
            %     sr = sr/max(sr); % normalize
            
            %%%% Detection Above Noise -Level
%             [nlvl, nstd] = noise_lvl_est(sr, prc); % cancelled for better time performance
            [nlvl, nstd] = deal(NoiseLevel);
            sum_pkTH = max(nlvl+2.5*nstd,  minTHLvl);
            sum_pkPR = max(2.5*nstd,  minTHLvl); % nlvl;
%             [pks_sum, ipks_sum ] = findpeaks(sr, 'MinPeakDistance',sum_pkDist,  'MinPeakHeight',sum_pkTH, 'MinPeakProminence', sum_pkPR);
            [pks_sum, ipks_sum ] = myFindPeaks(sr, sum_pkTH, sum_pkDist ); % Changed Sep2022
            %%%% OUTPUT
            FilterBank_Out.Active_Channels = active_channels;
            FilterBank_Out.Sum_out = sr;
            FilterBank_Out.det_pks = pks_sum;
            FilterBank_Out.det_ipks = ipks_sum;
            tsum = tout; % for the plot
        
        case "channel_detection"
            % Relevant Frequencies of the call 
            active_channels = ~isnan([FilterBank_Out.Ref_ipks]);
            nch = sum(active_channels);
            
            % all_detections
            all_ipks = [FilterBank_Out.Channels.ipks_shift];
            
            %%% histogram of detections
            % Gaussian Filter
            L =64;
            if ~clutterMode
                sigma = 30e-6*fs; %30e-6*fs; %7.5e-6*fs; % 7.5e-6*fs;  % 7.5e-6*fs; % 50e-6*fs; % the resolutin of detections 50 microsec %%% Change for Classifier
                sum_pkTH =  0.25*nch;   %  0.2*nch;% /3;
            else % 
                sigma = 5e-6*fs; %2.5e-6*fs; % 5e-6*fs;
                sum_pkTH =  max(0.1*nch, 1.5); %max(0.1*nch, 2.5)   %  0.2*nch;% /3;
            end % if ~clutterMode
            alpha = (L-1)/(2*sigma);
            g_win = gausswin(L,alpha);
                
            if ~isempty(all_ipks)
                bin_size = 1 ; %2 ; % 4; % 1; % 50e-6*fs;% 0.1/fs1e-3*fs;
                [N,edges] = histcounts(all_ipks,'BinWidth',bin_size);
                if bin_size == 1
                    mid_bins = edges(1:end-1);
                else
                    mid_bins = edges(1:end-1) + diff(edges)/2;
                end % if
                
                % Apply the filter
                smoothN = conv(N, g_win,"same");
                tsum = mid_bins/fs; % for the plot
                sigma = min(sigma, numel(smoothN)-2);
            else % if ~isempty(all_ipks)
                % if there ar no detections- zeros
                smoothN = zeros(size(FilterBank_Out.Channels(ifilt).response));
                tsum = ( 1:numel(smoothN) ) / fs; %sec;
                mid_bins = [];
                bin_size = 0;
            end % if ~isempty(all_ipks)

            %Final Detections    
            sum_pkDist = max(1.5*sigma, bin_size/2);
            if sum_pkDist >= numel(smoothN) -1
                sum_pkDist = numel(smoothN) -2;
            end
            sum_prTH = 0.25*sum_pkTH; % 0.1*sum_pkTH;


            if numel(smoothN) > 3
%                 [pks_sum, ipks_sum] = findpeaks(smoothN, 'MinPeakDistance', sigma,  'MinPeakHeight', sum_pkTH);
%                 [pks_sum, ipks_sum, ~, prs] = findpeaks(smoothN, 'MinPeakDistance', sum_pkDist,  'MinPeakHeight', sum_pkTH, 'MinPeakProminence', sum_prTH);
                [pks_sum, ipks_sum] = myFindPeaks(smoothN, sum_pkTH, sum_pkDist ); % Changed Sep2022
            elseif max(smoothN) > sum_pkTH 
                [pks_sum, ipks_sum] = max(smoothN);
            else
                pks_sum  = [];
                ipks_sum = [];
            end % if numel(smoothN) > 3
            
            % OnlyPositive ipks
            det_ipks = round(mid_bins(ipks_sum));

            %%%% OUTPUT
            FilterBank_Out.Active_Channels = active_channels;
            FilterBank_Out.Sum_out = sr;
            FilterBank_Out.detection_TH = sum_pkTH;
            FilterBank_Out.det_pks = pks_sum(det_ipks > 0);
            FilterBank_Out.det_ipks = det_ipks(det_ipks >0 ); % mid_bins(ipks_sum);
            FilterBank_Out.time_vec = tsum;
            FilterBank_Out.Sum_ipks_hist = smoothN;
    end % switch
    
    %% Figure 
    if debug_flag
        
        all_shifted_responses = reshape([FilterBank_Out.Channels.shifted_response], [], nch);
        hold(ax2,'off')
        plot(ax2, tout, all_shifted_responses, 'DisplayName', 'all channels')
        hold(ax2,'on')
        plot(ax2, [FilterBank_Out.Channels.ipks_shift]/fs, [FilterBank_Out.Channels.pks_shift],'.k', 'MarkerSize', 8,'DisplayName', 'pks')
        hold(ax3,'on')
        switch integration_method
            case  "summation"
                plot(ax3, tout, FilterBank_Out.Sum_out, 'k','LineWidth',2, 'DisplayName', 'Sum Shifted')
            case "channel_detection"
%                 plot(ax3, tout, FilterBank_Out.Sum_out, 'b','LineWidth',1, 'DisplayName', 'Sum Shifted')
%                 plot(ax3, tsum, FilterBank_Out.Sum_ipks_hist, 'k','LineWidth',2, 'DisplayName', 'Peaks histogram')
                try
                    bar(ax3, tsum, N, 'r', 'BarWidth', 1.5 , 'DisplayName', 'Peaks histogram')
                    plot(ax3, tsum, FilterBank_Out.Sum_ipks_hist, 'k','LineWidth',2, 'DisplayName', 'smooth Peaks histogram')
                catch
                    warning('No peaks Detectced')
                end
        end % switch
        plot(ax3,  FilterBank_Out.det_ipks/fs, FilterBank_Out.det_pks, 'og', 'DisplayName', 'Detections')
        plot(ax3, xlim, sum_pkTH*[1,1], 'r--')
        linkaxes([ax1, ax2, ax3],'x')
        xlim(ax2, [0 max(tout)]);
        title(ax2, 'Shifted channels')
        title(ax3, 'Summary')
        legend(ax3)
    end %
    
end % if ~reference_flag

end % main function

