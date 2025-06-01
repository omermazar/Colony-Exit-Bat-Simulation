function [resultsTable ] = ...
    main_jamMoth_test(test_params, AllParams, jamPreySignal, FB_params, target_params, plot_flags)

%% Input

%% Test_params - a user defindes sturct with hte field 'name' as name of the test

% %%%%% 'SingleEcho'
% Test single Echo with Jamming Noise with difdferns SIR
% test_params.name  = 'SingleEcho'; 
% test_params.sigType = "Search"; %  "Search" | "Approach" | "Buzz"
% test_params.sigLvl = 1; % [linear]
% test_params.SIR = [30 20 10 5 3 0 -3 -5 -10 -20]; % signal to Interference Ratio [dB]
% test_params.signalTime = [5 11 17 23 29 35]*1e-3; % startTime of the signal  to the jamming signal 30 ms to 
% test_params.jamSatrtTime = 10*1e-3;
% test_params.integration_method = "channel_detection"; %  { "summation", "channel_detection"}channel_detection

if ~isfield(test_params, 'detection_mode')
    test_params.detection_mode = 'actual';
end

% default params

%%  Filter Bank params

if nargin < 4
    FB_params.num_of_filters = 120; %120; %80;
    FB_params.fs = 200e3;
    FB_params.f_low = 20e3;
    FB_params.f_high = 80e3;
    FB_params.response_time = 0.5e-3; % 0.5e-3; % good for f_low = 20kHz
    FB_params.fft_resolution = 1e3;
    FB_params.f_lowpass= 12e3; % 8e3;
end % if nargin < 3

names = fieldnames(FB_params);
for i=1:length(names)
    eval([names{i} '=FB_params.' names{i} ]);
end

%%% Target Reference
if nargin < 6
    target_params.n_echoes = 2; 
    target_params.time_diff     = 0.08e-3; % 0.15e-3; % time differnce between the echoes of the target
    target_params.echo_lvl      = [1 1]; % the level of each echo of the target
    target_params.start_period  = 1e-3; % period of quiet before the signal (s)
    target_params.response_time = 0.6e-3; % the total time in th PS space (s)

end

names = fieldnames(target_params);
for i=1:length(names)
    eval([names{i} '=target_params.' names{i} ]);
end


%%% plot flags
if nargin < 7
    plot_flags.plot_tar_ref_flg  = 1;
    plot_flags.plot_norm_flag    = 1; % 1;
    plot_flags.plot_all_reponses = 1;
end

%%  Warnings
warning('off', 'signal:findpeaks:largeMinPeakHeight')
warning('off','MATLAB:legend:IgnoringExtraEntries')
% warning('off', 'Invalid MinPeakHeight')


%%%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
resultsTable = []; 

%% Input Signal (the bat call)
%%%%%  Chirp Params %%%%
switch test_params.sigType
    case 'Search'
        % Params for Seach call (BBB)
        Chirp_params.call_dur    = 5e-3; % 5e-3; %
        Chirp_params.PulsePower  = 1;
        Chirp_params.fchirp_high = 52e3; % 60e3;
        Chirp_params.fchirp_low  = 32e3;

    case 'Approach'
        % Params for Approach call (BBB)
        Chirp_params.call_dur    = 4e-3; % 5e-3; %
        Chirp_params.PulsePower  = 1;
        Chirp_params.fchirp_high = 78e3; % 60e3;
        Chirp_params.fchirp_low  = 32e3;

    case 'Buzz'
        % Params for Buzz call (BBB)
        Chirp_params.call_dur    = 1.5e-3; % 5e-3; %
        Chirp_params.PulsePower  = 1;
        Chirp_params.fchirp_high = 35e3; % 60e3;
        Chirp_params.fchirp_low  = 25e3;
end % switch test_params.sigType

names = fieldnames(Chirp_params);
for i=1:length(names)
    eval([names{i} '= Chirp_params.' names{i} ]);
end

%% %%% The Jamming signal
if nargin < 3
    jamPreySignal = generateAcousicPreyJamming(AllParams.PreyFlightParams, fs);
    jamPreySignal = jamPreySignal ./ 10.^(AllParams.PreyFlightParams.JamTxLvl/20);
%     jamPreySignal = jamPreySignal ./ std(jamPreySignal);
end 

%% Calculate the FilterBank Impulse Response
% LPF 
[LPF_filt.b_filt, LPF_filt.a_filt] = butter(6, f_lowpass/(fs/2)); % butter(6, f_lowpass/(fs/2));

[AllFiltersParams]  = FilterBank_Acous_Resp(num_of_filters, fs, f_low, f_high, response_time, fft_resolution);
fc  = AllFiltersParams.filter_fc;


%% Reference Input Signal (the bat call)
tin = [0:1/fs:call_dur]; % time vector
sig_in = PulsePower * chirp(tin, fchirp_high, max(tin), fchirp_low, 'logarithmic'); % chirp(tin, f_high, max(tin), f_low, 'logarithmic')

%%%%%% THE REFERENCE RESPONSE %%%%%
[FilterBank_Ref] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_in, fs,1, [], test_params.integration_method, 0);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The Test
% sum_tbl = struct2table(target_results); 

eco_time = 10*1e-3*fs; % the start time of the eco

switch test_params.name
    case {'SingleEcho'}
        %%
        % Build the input signal 
        techo    = test_params.signalTime(1); % the time of first echo
        tend     = 2e-3;  % quiet time after the end of the last echo (s)
        snrdB    = 50; % Signal to Noise Ratio   
        echo_lvl = test_params.sigLvl ./ std(sig_in);
        tot_time = 50e-3; % the total time of the detections
        
        noiselvl = max(echo_lvl) / 10.^(snrdB/20);
        t = (0:1/fs:tot_time)*1e3; % ms
        sig_echo = zeros(size(t));
        sig_jam = noiselvl * randn(size(t));
        nEchoes = numel(test_params.signalTime);
        time_diff = diff(test_params.signalTime);
        
        detectionTimeTol = 50e-6;

         %% The Response without interfewrence
        sig_echo1 = multi_echo_gen(sig_in, nEchoes, time_diff, echo_lvl, snrdB, fs, techo, tend);
        sig_echo(1:numel(sig_echo1)) = sig_echo1;

        FilterBank_Echoes = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_echo, fs, false, FilterBank_Ref, test_params.integration_method, 1);
        sgtitle('Reference')
        detected_echoesOnly = FilterBank_Echoes.det_ipks/ fs;
        
        ref_echoesTimes = test_params.signalTime;
%         [k,dist] = dsearchn(ref_echoesTimes', detected_echoesOnly');
%         ixDet = dist <= detectionTimeTol;
        results.Reference = my_detection_compare(ref_echoesTimes, FilterBank_Echoes, detectionTimeTol, fs);
       
        %  %% Output
        callType = test_params.sigType;
        refTable = restults2Table(results.Reference, "Reference", nEchoes, inf, callType);
        resultsTable = [resultsTable; refTable]; 
        
        %% The Jamming Response (inreference Only)    
        jamStart = round(test_params.jamSatrtTime * fs);
        idx_jam = jamStart + (1:numel(jamPreySignal) );
        sig_jam(idx_jam) = jamPreySignal;
        
        [FilterBank_jam] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_jam, fs, false, FilterBank_Ref, test_params.integration_method, 1);
        results.JammingOnly = my_detection_compare(ref_echoesTimes, FilterBank_jam, detectionTimeTol, fs);
        sgtitle('Jamming')

        %  %% Output
        tempTable = restults2Table(results.JammingOnly, "JamOnly", nEchoes, -inf, callType);
        resultsTable = [resultsTable; tempTable]; 

        %% The Rx Signal (echoes + Jaming)
        for SIR = test_params.SIR
            sig_jam = noiselvl * randn(size(t));
            sig_jam(idx_jam) = jamPreySignal ./ 10.^(SIR/20);
            sig_tot = sig_echo + sig_jam;
            [FilterBank_tot] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_tot, fs, false, FilterBank_Ref, test_params.integration_method, 0);
            results.Total = my_detection_compare(ref_echoesTimes, FilterBank_tot, detectionTimeTol, fs);

            % figure
            fig = plotJamTest(t, sig_echo, sig_jam, sig_tot, FilterBank_tot, fs, test_params, results.Total, ref_echoesTimes);
            sgtitle(fig, [test_params.sigType, ', SIR: ', num2str(SIR), ' dB'])

            %summary   %% Output
            tempName = strcat("SIR", num2str(SIR), "dB");
            tempTable = restults2Table(results.Total, tempName, nEchoes, SIR, callType);
            resultsTable = [resultsTable; tempTable]; 

             save_flag = 1;
             if save_flag
                 file_name_prefix = strcat('\', test_params.name, '_', test_params.integration_method, '_', test_params.sigType, '_') ;
                 save_path = 'D:\Dropbox\University\מחקר\experiments\Moth_Jamming\Figures';

                 %%%% figures
                 filename = strcat(save_path, file_name_prefix, 'SNR_', num2str(SIR));
                 savefig (fig,  filename);
                 print(fig, filename, '-djpeg', '-r0', '-noui')   
             end % if save_flag
        end % for SIR = test_params.SIR

        %% end SingleEcho

      
        
    end % switch


%% summary
%%%% Add Chirp Params to sum_table
%% save_data
save_flag = 1;
if save_flag
   save_path = 'D:\Dropbox\University\מחקר\experiments\Moth_Jamming\Results';
   file_name_prefix = strcat('\', test_params.name, '_', test_params.integration_method, '_', test_params.sigType, '_') ;
   save( strcat(save_path, file_name_prefix, 'sumResults'), 'resultsTable')
end % if save_flag
