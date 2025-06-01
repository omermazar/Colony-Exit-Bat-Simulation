function [sum_tbl, sum_pks, FR_all_tests, FilterBank_target, sig_target] = main_targetClassifierTest(test_params, FB_params, Chirp_params, target_params, plot_flags)

%% Input

% Test_params - a user defindes sturct with hte field 'name' as name of the test
%%%%% 'delay_sensitivity'
% test_params.name = 'delay_sensitivity'
% test_params.delays_ms = [-0.1e-3: 5/fs: 0.15e-3 ];
% test_params.SNR = 30; % dB
% 
% %%%%% 'SNR"
% test_params.name = 'SNR'
% test_params.SNR = [40, 30, 20, 10, 6 , 3 , 0 , -3 ]; % dB
% test_params.detection_mode = 'actual'; % 'actual' ; 'FB_detection'

%%%%%% 'OnePoint_target'
% test_params.name = OnePoint_target'; %'OnePoint_target_new'
% test_params.SNR = [40 30 20 15 10 6 3 0 -3 -6 -10 -20]; % dB

%%%%%% 'Dipoles'
% test_params.name = 'Dipoles'
% test_params.dipole_diffs = [1 0.6 0.4 0.3 0.2 0.175 0.16 0.14 0.13 0.12 0.1 0.08 0.05]*1e-3; % sec
% test_params.SNR = 30; % dB;
% test_params.detection_mode = 'actual'

%%%%%% 'Clutter'
% test_params.name = 'Clutter';
% test_params.target2clutter_time_diff = 1e-3; %sec;
% test_params.echoes_num = 10; % N
% test_params.mean_diff = 0.12e-3; % sec;
% test_params.std_diff = 0.06e-3; % sec;
% test_params.mean_lvl = 1; % no-units
% test_params.std_lvl = 0.2; % 
% test_params.SNR = 30;
% test_params.tests_num = 1;

%%%%% 'Interference'
% test_params.name = 'Interference';
% test_params.freq_diff = [-2000 -1000 -500 -250 100 0 3000 500 1000 3000 5000]; % [-1, -0.5, 0.5, 1, 3, 5]*1e3; 
% test_params.bw_diff = 0; % sec; 
% test_params.par2tests = 'freq_diff' ; % 'freq_diff' , 'bw_diff'
% test_params.SNR = 30;

%%%%% 'DOA'
% test_params.name = 'DOA'; 
% test_params.angle2target = [0 5 10 20 30 45 60 90] ; % degrees
% test_params.SNR = 30; % dB

%%%%% 'random_noise'
% test_params.name = 'random_noise'; 
% test_params.noise_lvl = 0 ; %dB
% test_params.SNR = 30; % dB
% test_params.tests_num = 10;

%%%%%% 'ChirpBw'
% test_params.name = 'ChirpBw'
% test_params.chirpBW = [5 7 8 10 15 20 30 40]*1e3;
% test_params.SNR = 30; % dB
% test_params.detection_mode = 'actual'

if ~isfield(test_params, 'detection_mode')
    test_params.detection_mode = 'actual';
end

% default params

%%%%%  Filter Bank params

if nargin < 2
    FB_params.num_of_filters = 120; %120; %80;
    FB_params.fs = 200e3;
    FB_params.f_low = 20e3;
    FB_params.f_high = 80e3;
    FB_params.response_time = 0.5e-3; % 0.5e-3; % good for f_low = 20kHz
    FB_params.fft_resolution = 1e3;
    FB_params.f_lowpass= 12e3; % 8e3;
end % if nargin < 1

names = fieldnames(FB_params);
for i=1:length(names)
    eval([names{i} '=FB_params.' names{i} ]);
end

%%%%%  Chirp Params %%%%
if nargin < 3
    Chirp_params.call_dur = 7e-3; %4e-3; % 5e-3; %
    Chirp_params.PulsePower = 1;
    Chirp_params.fchirp_high = 50e3; % 60e3;
    Chirp_params.fchirp_low = 32e3; %40e3;
end % if nargin < 2

names = fieldnames(Chirp_params);
for i=1:length(names)
    eval([names{i} '= Chirp_params.' names{i} ]);
end

%%% Target Reference
if nargin < 4
    target_params.n_echoes      = 2; 
    target_params.time_diff     = 0.08e-3; % 0.15e-3; % time differnce between the echoes of the target
    target_params.echo_lvl      = [1 1]; % the level of each echo of the target
    target_params.start_period  = 1e-3; % period of quiet before the signal (s)
    target_params.response_time = 0.6e-3; % the total time in th PS space (s)

end

names = fieldnames(target_params);
for i=1:length(names)
    eval([names{i} '=target_params.' names{i} ]);
end

%%% plot flaggs
if nargin < 5
    plot_flags.plot_tar_ref_flg  = 1;
    plot_flags.plot_norm_flag    = 1; % 1;
    plot_flags.plot_all_reponses = 1;
end

%%  Warnings
warning('off', 'signal:findpeaks:largeMinPeakHeight')
warning('off','MATLAB:legend:IgnoringExtraEntries')
% warning('off', 'Invalid MinPeakHeight')

%% Calculathe the FilterBank Impulse Response
% LPF 
[LPF_filt.b_filt, LPF_filt.a_filt] = butter(6, f_lowpass/(fs/2)); % butter(6, f_lowpass/(fs/2));

[AllFiltersParams]  = FilterBank_Acous_Resp(num_of_filters, fs, f_low, f_high, response_time, fft_resolution);
fc = AllFiltersParams.filter_fc;


%% Reference Input Signal
tin = [0:1/fs:call_dur]; % time vector
sig_in = PulsePower * chirp(tin, fchirp_high, max(tin), fchirp_low, 'logarithmic'); % chirp(tin, f_high, max(tin), f_low, 'logarithmic')

%%%%%% THE REFERENCE RESPONSE %%%%%
[FilterBank_Ref] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_in, fs,1);

%% Reference Target Signature
n_echoes =  target_params.n_echoes;
echo_lvl = target_params.echo_lvl;
snrdB    = 100;
tstart   = start_period;
tend     = start_period;
sig_target = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, snrdB, fs, tstart, tend);


%% The Target Response
[FilterBank_target] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_target, fs, 0, FilterBank_Ref);
 
active_ind = FilterBank_target.Active_Channels;
n_active = sum(active_ind);

nresponse = target_params.response_time * fs;
nstart = tstart * fs;

mat_target_full =  reshape([FilterBank_target.Channels(active_ind).shifted_response], [], n_active);
ix_mat_target = round((nstart-nresponse):(nstart+nresponse));
mat_target = mat_target_full( ix_mat_target, : );
ix_start_target = round(nresponse+1);
FR_target = mat_target(ix_start_target, : );
f = fc(active_ind);

% for reference
[target_results, pks_target, time_pks_target] = analyze_echo2ref(sig_target, sig_target, mat_target, FR_target, mat_target, FR_target, true, [], [], f, fs);
target_results.name            = "Reference";
target_results.changed_param   = nan;
target_results.delay_error     = 0; 
target_results.is_detected_flg = 0;
target_results.snrdB           = snrdB;
pks_target.name = "Reference";
pks_target.changed_param = nan;


if plot_flags.plot_tar_ref_flg
    my_plot_tar_ref(mat_target, FR_target, ix_start_target, f, fs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The Test
sum_tbl = struct2table(target_results); 


sum_pks = pks_target; 
eco_time = 10*1e-3*fs; % the start time of the eco

switch test_params.name
    case {'SNR', 'SNR_new'}
        %%
        detection_num = numel(test_params.SNR);
        FR_all_tests = nan(detection_num, numel(f));
        irow = 1;
        
        % the echo parameters
        n_echoes = target_params.n_echoes;
        echo_lvl = target_params.echo_lvl;
        snrdB    = nan; % changed each test
        techo    = 10e-3;
        tend     = 2e-3;
        
        nstart = round(techo * fs);
        ix_mat_echo = round((nstart-nresponse):(nstart+nresponse));  % the index of the echo respons
        
        for SNR = test_params.SNR
            snrdB = SNR;
            % the echo
            sig_echo = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, snrdB, fs, techo, tend);
            % The response
            [FilterBank_echo] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_echo, fs, false, FilterBank_Ref);
            mat_full = reshape([FilterBank_echo.Channels(active_ind).shifted_response], [], n_active);          
            
            % take the relevant echo (from detection)
            actual_echo_time = round((techo - tstart)*fs)+1;
            switch test_params.detection_mode
                case 'actual'
                    det_time = actual_echo_time;
                    ndelay = 0;
                    is_detected_flg = 1;
                case 'FB_detection'
                    if any(FilterBank_echo.det_ipks)
                        [k, fb_det] = dsearchn(FilterBank_echo.det_ipks'/fs, (techo - tstart) );
                        det_time = FilterBank_echo.det_ipks(k) - tstart*fs +1;
                        ndelay = det_time - actual_echo_time;
                        is_detected_flg = 1;
                    else
                        det_time = nan;
                        ndelay = det_time - actual_echo_time;
                        is_detected_flg = 0;
                    end % if any(FilterBank_echo.det_ipks)


            end % switch test_params.detection_mode
            % the error of detection
            delay_error = ndelay / fs;
            
             display(['SNR: ' , num2str(SNR), ' times: ', num2str(FilterBank_echo.det_ipks/fs,4) , ' delay_error: ', num2str(delay_error,4)])
                 
             if is_detected_flg
                 % the detected signals (shifted)
                 mat_echo = mat_full( (ix_mat_echo + ndelay), : );
                 FR_echo = mat_echo(round(nresponse+1), : );
                 sig_echo = sig_echo( det_time:(det_time + numel(sig_target)) );
                 % the result
                 [test_results, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                     mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = analyze_echo2ref(sig_target, sig_echo, mat_target, FR_target, mat_echo, FR_echo, ...
                     false, pks_target, time_pks_target, f, fs );

             else % (not detected) if is_detected_flg
                 % fill nan values
                 sig_echo = sig_echo( actual_echo_time : (actual_echo_time + numel(sig_target)) ); % just for plot
                 [mat_echo, FR_echo, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                     mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2] = deal(nan);
                 test_results = fill_struct_fields(target_results, nan, {'is_detected_flg' , 'name' , 'changed_param', 'delay_error'});
                 pks_out = fill_struct_fields(pks_target, nan, { 'name' , 'changed_param',});
                 time_pks_out = fill_struct_fields(time_pks_target, nan, { 'name' , 'changed_param',});

             end % if is_detected_flg
            
            test_results.is_detected_flg = is_detected_flg;
            test_results.delay_error = delay_error;
            test_results.name = string(test_params.name);
            test_results.changed_param = SNR;
            pks_out.name = string(test_params.name);
            pks_out.changed_param = SNR;
            
            %%%% plot %%%%%
            if plot_flags.plot_norm_flag & is_detected_flg
                my_plot_norm_results(test_results,  sig_target, sig_echo, pks_out, pks_target, time_pks_out, time_pks_target, mat_target, mat_echo, FR_target, FR_echo,   ...
                    mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                    mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2, f, fs);
            end % plot

            %%% update summary table with carreun results
            sum_tbl = [sum_tbl; struct2table(test_results)];
            sum_pks = [sum_pks; pks_out];
            %%% update all Frequecy Reponses
            FR_all_tests(irow,:) = FR_echo;
            irow = irow+1;

        end % for SNR
        
        sum_pks = struct2table(sum_pks);
%         %%%% plot %%%%%
%         if plot_flags.plot_all_reponses
%             my_plot_all_freponses(test_params, FR_target, FR_all_tests, f, fs);
%         end % plot

        %% end SNR Test
    
    case {'delay_sensitivity', 'delay_sensitivity30dB'}
        %%
        snrdB = test_params.SNR;
        techo = tstart;
        sig_echo = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, snrdB, fs, techo, tend); % same as target except SNR
        
%        FilterBank_echo = FilterBank_target;
        [FilterBank_echo] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_echo, fs, false, FilterBank_Ref);
        mat_echo_full =  reshape([FilterBank_echo.Channels(active_ind).shifted_response], [], n_active);
        actual_echo_time = round((techo - tstart)*fs)+1;
        is_detected_flg = 1; % in this test there are no detections

        % the tests
        ndelays =  round(test_params.delays_ms * fs);
        FR_all_tests = nan(numel(ndelays), numel(f));
        irow =1;
        for idelay = ndelays

            mat_echo = mat_echo_full(ix_mat_target + idelay, :);
            FR_echo = mat_echo(ix_start_target + idelay, : );
            
            [test_results, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = analyze_echo2ref(sig_target, sig_echo,mat_target, FR_target, mat_echo, FR_echo, ...
                false, pks_target, time_pks_target, f, fs );
            
            test_results.is_detected_flg = is_detected_flg;
            test_results.delay_error = idelay/fs;
            test_results.name = string(test_params.name);
            test_results.changed_param = idelay*1000/fs;
            pks_out.name = string(test_params.name);
            pks_out.changed_param = idelay*1000/fs;

            %%%% plot %%%%%
            if plot_flags.plot_norm_flag
               my_plot_norm_results(test_results,  sig_target, sig_echo, pks_out, pks_target, time_pks_out, time_pks_target, mat_target, mat_echo, FR_target, FR_echo,   ...
                    mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                    mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2, f, fs);
            end % plot

            %%% update summary table with carreun results
            sum_tbl = [sum_tbl; struct2table(test_results)];
            sum_pks = [sum_pks; pks_out];
            %%% update all Frequecy Reponses
            FR_all_tests(irow,:) = FR_echo;
            irow = irow+1;
        
        end % for idelay = ndelays
        
        sum_pks = struct2table(sum_pks);

        %%%% plot %%%%%
%         if plot_flags.plot_all_reponses
%             my_plot_all_freponses(test_params, FR_target, FR_all_tests, f, fs);
%         end % plot
        %% end of 'delay_sensitivity'

    case {'OnePoint_target', 'OnePoint_target_new'}
         %%
        detection_num = numel(test_params.SNR);
        FR_all_tests = nan(detection_num, numel(f));
        irow = 1;
        
        % the echo parameters
        n_echoes = 1;
        echo_lvl = 1;
        snrdB    = nan; % changed each test
        techo   = 10e-3;
        tend     = 2e-3;
        time_diff = 0; 

        nstart = round(techo * fs);
        ix_mat_echo = round((nstart-nresponse):(nstart+nresponse));  % the index of the echo respons
        
        for SNR = test_params.SNR
            snrdB = SNR;
            % the echo
            sig_echo = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, snrdB, fs, techo, tend);
            % The response
            [FilterBank_echo] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_echo, fs, false, FilterBank_Ref);
            
            display(['SNR: ' , num2str(SNR), ' times: ', num2str(FilterBank_echo.det_ipks/fs,4) ])

            mat_full = reshape([FilterBank_echo.Channels(active_ind).shifted_response], [], n_active);          
            mat_echo = mat_full( ix_mat_echo, : );
            FR_echo = mat_echo(round(nresponse+1), : );
             % take the relevant echo (from detection)
            det_time =  round((techo - tstart)*fs)+1;
            sig_echo = sig_echo( det_time:(det_time + numel(sig_target)) );

            [test_results, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = analyze_echo2ref(sig_target, sig_echo, mat_target, FR_target, mat_echo, FR_echo, ...
                false, pks_target, time_pks_target, f, fs );
            
            test_results.name = string(test_params.name);
            test_results.changed_param = SNR;
            pks_out.name = string(test_params.name);
            pks_out.changed_param = SNR;
            
            %%%% plot %%%%%
            if plot_flags.plot_norm_flag
               my_plot_norm_results(test_results,  sig_target, sig_echo, pks_out, pks_target, time_pks_out, time_pks_target, mat_target, mat_echo, FR_target, FR_echo,   ...
                    mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                    mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2, f, fs);
            end % plot

            %%% update summary table with carreun results
            sum_tbl = [sum_tbl; struct2table(test_results)];
            sum_pks = [sum_pks; pks_out];
            %%% update all Frequecy Reponses
            FR_all_tests(irow,:) = FR_echo;
            irow = irow+1;

        end % for OnePoint

        sum_pks = struct2table(sum_pks);

        

        %% end OnePoint_target
        case {'Dipoles', 'Dipoles_vs_Echoes'}
        %%
        detection_num = numel(test_params.SNR);
        FR_all_tests = nan(detection_num, numel(f));
        irow = 1;
        
        % the echo parameters
        n_echoes  = 2;
        echo_lvl  = 1;
        snrdB     = test_params.SNR;
        techo     = 10e-3;
        tend      = 2e-3;
        time_diff = nan; % chaneged each test 

        nstart = round(techo * fs);
        ix_mat_echo = round((nstart-nresponse):(nstart+nresponse));  % the index of the echo respons
        
        for time_diff = test_params.dipole_diffs
            % the echo
            sig_echo = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, snrdB, fs, techo, tend);
             actual_echo_time = round((techo - tstart)*fs)+1;
            % The response
            [FilterBank_echo] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_echo, fs, false, FilterBank_Ref);
            mat_full = reshape([FilterBank_echo.Channels(active_ind).shifted_response], [], n_active);  

            % take the relevant echo (from detection)
            switch test_params.detection_mode
                case 'actual'
                    det_time = actual_echo_time;
                    ndelay = 0;
                    is_detected_flg = 1;
                case 'FB_detection'
                    if any(FilterBank_echo.det_ipks)
                        [k, fb_det] = dsearchn(FilterBank_echo.det_ipks'/fs, (techo - tstart) );
                        det_time = FilterBank_echo.det_ipks(k) - tstart*fs +1;
                        ndelay = det_time - actual_echo_time;
                        is_detected_flg = 1;
                    else
                        det_time = nan;
                        ndelay = det_time - actual_echo_time;
                        is_detected_flg = 0;
                    end % if any(FilterBank_echo.det_ipks)

            end % switch test_params.detection_mode
             % the error of detection
            delay_error = ndelay / fs;

            display(['Dipole: ' , num2str(time_diff), ' times: ', num2str(FilterBank_echo.det_ipks/fs,4) , ' delay_error: ', num2str(delay_error,4)])

%             
%             mat_echo = mat_full( ix_mat_echo, : );
%             FR_echo = mat_echo(round(nresponse+1), : );
%              % take the relevant echo (from detection)
%             det_time =  round((techo - tstart)*fs)+1;
%             sig_echo = sig_echo( det_time:(det_time + numel(sig_target)) );
            if is_detected_flg
                 % the detected signals (shifted)
                 mat_echo = mat_full( (ix_mat_echo + ndelay), : );
                 FR_echo = mat_echo(round(nresponse+1), : );
                 sig_echo = sig_echo( det_time:(det_time + numel(sig_target)) );
                 % the result
                 [test_results, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                     mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = analyze_echo2ref(sig_target, sig_echo, mat_target, FR_target, mat_echo, FR_echo, ...
                     false, pks_target, time_pks_target, f, fs );

             else % (not detected) if is_detected_flg
                 % fill nan values
                 sig_echo = sig_echo( actual_echo_time : (actual_echo_time + numel(sig_target)) ); % just for plot
                 [mat_echo, FR_echo, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                     mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2] = deal(nan);
                 test_results = fill_struct_fields(target_results, nan, {'is_detected_flg' , 'name' , 'changed_param', 'delay_error'});
                 pks_out = fill_struct_fields(pks_target, nan, { 'name' , 'changed_param',});
                 time_pks_out = fill_struct_fields(time_pks_target, nan, { 'name' , 'changed_param',});

             end % if is_detected_flg

%             [test_results, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
%                 mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = analyze_echo2ref(sig_target, sig_echo, mat_target, FR_target, mat_echo, FR_echo, ...
%                 false, pks_target, time_pks_target, f, fs );
            
            test_results.is_detected_flg = is_detected_flg;
            test_results.delay_error = delay_error;
            test_results.name = string(test_params.name);
            test_results.changed_param = time_diff;
            pks_out.name = string(test_params.name);
            pks_out.changed_param = time_diff;
            
            %%%% plot %%%%%
            if plot_flags.plot_norm_flag
                my_plot_norm_results(test_results,  sig_target, sig_echo, pks_out, pks_target, time_pks_out, time_pks_target, mat_target, mat_echo, FR_target, FR_echo,   ...
                    mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                    mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2, f, fs);
            end % plot

            %%% update summary table with carreun results
            try
                sum_tbl = [sum_tbl; struct2table(test_results)];
                sum_pks = [sum_pks; pks_out];
            catch
                warning(['Table Conc Error time_diff = ', num2str(time_diff)]);
            end % try 
            %%% update all Frequecy Reponses
            FR_all_tests(irow,:) = FR_echo;
            irow = irow+1;

        end % for time_diff
        
        sum_pks = struct2table(sum_pks);

        %% end Dipoles
    case 'Clutter'
        %%
        % add clutter data to summary table
        sum_tbl.clutter_diffs = nan(1,test_params.echoes_num );
        sum_tbl.clutter_lvls = nan(1,test_params.echoes_num );

        for iTest = 1: test_params.tests_num % number of repetions
            % create the clutter
            ix_negative = true(1,test_params.echoes_num);
            clutter.diff = nan(1,test_params.echoes_num);
            clutter.lvl = abs( test_params.mean_lvl + test_params.std_lvl*randn(1, test_params.echoes_num) );
            while any(ix_negative)
                clutter.diff(ix_negative) = test_params.mean_diff + test_params.std_diff*randn(1, sum(ix_negative));
                ix_negative = clutter.diff <= 0;
            end
            clutter.diff(1) =  0; % no meaning fo the first

            % create the echo from the clutter
            % the echo parameters
            n_echoes  = test_params.echoes_num;
            echo_lvl  = clutter.lvl;
            snrdB     = test_params.SNR;
            techo     = 10e-3 + test_params.target2clutter_time_diff;
            tend      = 5e-3;
            time_diff = clutter.diff(2:end);
            
            clutter_echo = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, snrdB, fs, techo, tend);

            %%% The FilterBank Repsonse to Clutter
            [FilterBank_clutter] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, clutter_echo, fs, false, FilterBank_Ref);

            %%% look at tre resuts of the detected lutter
            detect_ntimes = FilterBank_clutter.det_ipks;
            detection_num = numel(detect_ntimes);
            FR_all_tests = nan(detection_num, numel(f));
            irow = 1;

            nstart = round(techo * fs);
            %         ix_mat_echo = round((nstart-nresponse):(nstart+nresponse));  % the index of the echo respons
            %

            for k = 1:detection_num
                % the echo
                %             sig_echo = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, snrdB, fs, techo, tend);
                %             % The response
                %             [FilterBank_echo] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_echo, fs, false, FilterBank_Ref);

                mat_full = reshape([FilterBank_clutter.Channels(active_ind).shifted_response], [], n_active);
                % the relvant index of each dectcion
                ix_mat_echo = round((detect_ntimes(k)-nresponse):(detect_ntimes(k)+nresponse));

                mat_echo = mat_full( ix_mat_echo, : );
                FR_echo = mat_echo(round(nresponse+1), : );
                 % take the relevant echo (from detection)
                det_time =  round((techo - tstart)*fs)+1;
%                   sig_echo = sig_echo( det_time:(det_time + numel(sig_target)) );
                sig_echo = clutter_echo( det_time:(det_time + numel(sig_target)) );

                [test_results, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                    mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = analyze_echo2ref(sig_target, sig_echo, mat_target, FR_target, mat_echo, FR_echo, ...
                    false, pks_target, time_pks_target, f, fs );

                test_results.name             = string(test_params.name);
                test_results.changed_param    = k;
                test_results.is_detected_flg  = nan;
                test_results.delay_error      = nan; 
                test_results.clutter_diffs    = clutter.diff;
                test_results.clutter_lvls     = clutter.lvl;
                pks_out.name                  = string(test_params.name);
                pks_out.changed_param         = k;

                %%%% plot %%%%%
%                 sig_echo = clutter_echo;
                if plot_flags.plot_norm_flag
                    my_plot_norm_results(test_results,  sig_target, sig_echo, pks_out, pks_target, time_pks_out, time_pks_target, mat_target, mat_echo, FR_target, FR_echo,   ...
                    mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                    mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2, f, fs);
                end % plot

                %%% update summary table with carreun results
                sum_tbl = [sum_tbl; struct2table(test_results)];
                sum_pks = [sum_pks; pks_out];
                %%% update all Frequecy Reponses
                FR_all_tests(irow,:) = FR_echo;
                irow = irow+1;

            end % for k

        end % for iTest = 1: test_params.tests_num

        sum_pks = struct2table(sum_pks);
        %% end 'Clutter'

    case 'Interference'
        %%
        tested_parameter = test_params.(test_params.par2tests);
        FR_all_tests = nan(numel(tested_parameter), numel(f));
        irow = 1;
        
        % the echo parameters
        n_echoes = 2;
        echo_lvl = 1;
        snrdB    = test_params.SNR; 
        techo   = 10e-3;
        tend     = 2e-3;
        
        % the interference signal parameter
        bw = Chirp_params.fchirp_high - Chirp_params.fchirp_low + test_params.bw_diff;
        Int_params.call_dur = Chirp_params.call_dur; % 5e-3; %
        Int_params.PulsePower = Chirp_params.PulsePower;
        Int_params.fchirp_high = Chirp_params.fchirp_low + test_params.freq_diff + bw; % 60e3;
        Int_params.fchirp_low =  Chirp_params.fchirp_low + test_params.freq_diff; 
        if numel( Int_params.fchirp_low) == 1
             Int_params.fchirp_low =  Int_params.fchirp_low*ones(size( Int_params.fchirp_high));
        end % if numel()
        nstart = round(techo * fs);
        ix_mat_echo = round((nstart-nresponse):(nstart+nresponse));  % the index of the echo respons
        
        for k = 1:numel(tested_parameter)
            % the interferference
            tint = [0:1/fs:Int_params.call_dur]; % time vector
            int_in = PulsePower * chirp(tint, Int_params.fchirp_high(k), max(tint), Int_params.fchirp_low(k), 'logarithmic'); 
           
            % the echo refelcted by target from interfernce
            sig_echo = multi_echo_gen(int_in, n_echoes, time_diff, echo_lvl, snrdB, fs, techo, tend);

            % The response
            [FilterBank_echo] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_echo, fs, false, FilterBank_Ref);
            display([test_params.par2tests, ': ' , num2str(tested_parameter(k)), ' times: ', num2str(FilterBank_echo.det_ipks/fs,4) ])
            mat_full = reshape([FilterBank_echo.Channels(active_ind).shifted_response], [], n_active);          
            mat_echo = mat_full( ix_mat_echo, : );
            FR_echo = mat_echo(round(nresponse+1), : );
            % take the relevant echo (from detection)
            det_time =  round((techo - tstart)*fs)+1;
            sig_echo = sig_echo( det_time:(det_time + numel(sig_target)) );

            [test_results, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = analyze_echo2ref(sig_target, sig_echo, mat_target, FR_target, mat_echo, FR_echo, ...
                false, pks_target, time_pks_target, f, fs );
            
            test_results.name = string(test_params.name);
            test_results.changed_param = tested_parameter(k);
            pks_out.name = string(test_params.name);
            pks_out.changed_param = tested_parameter(k);
            
            %%%% plot %%%%%
            if plot_flags.plot_norm_flag
                my_plot_norm_results(test_results, sig_target, sig_echo, pks_out, pks_target, time_pks_out, time_pks_target, mat_target, mat_echo, FR_target, FR_echo,   ...
                    mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                    mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2, f, fs);
            end % plot

            %%% update summary table with carreun results
            sum_tbl = [sum_tbl; struct2table(test_results)];
            sum_pks = [sum_pks; pks_out];
            %%% update all Frequecy Reponses
            FR_all_tests(irow,:) = FR_echo;
            irow = irow+1;

        end % for SNR
        
        sum_pks = struct2table(sum_pks);
%         %%%% plot %%%%%
%         if plot_flags.plot_all_reponses
%             my_plot_all_freponses(test_params, FR_target, FR_all_tests, f, fs);
%         end % plot

        %% end Interference Test
        
        case 'DOA'
        %%
        detection_num = numel(test_params.angle2target);
        FR_all_tests = nan(detection_num, numel(f));
        irow = 1;
        
        % the echo parameters
        n_echoes = 2;
        echo_lvl = 1;
        snrdB    = test_params.SNR; 
        techo    = 10e-3;
        tend     = 2e-3;
        
        
         % the freqs of the echo (the transmitted signal)
         bw = Chirp_params.fchirp_high - Chirp_params.fchirp_low;
         A = fchirp_high;
         b = log10(fchirp_low/fchirp_high)./numel(sig_in);
         FreqsVec = A*10.^(b*(1:numel(tin)) );
         % % %        [FreqsVec, ~] = CalcFreqsForPulse(Chirp_params.call_dur, Chirp_params.fchirp_low, bw, 'Logarithmic' )
         
        % paramters for Directivity calculation
        % CalculateSignalAtten(0.1, 0, 0.02, 20, 'meters' );
        Distances = 0.5;
        freq = FreqsVec/1000; % kHz
        CoordinationFlag =  'meters';
        TargetArea = 0.02^2*pi;

       
        nstart      = round(techo * fs);
        ix_mat_echo = round((nstart-nresponse):(nstart+nresponse));  % the index of the echo respons
        
        for teta = deg2rad(test_params.angle2target)
            % the diretctivity from each direction of each frequency
            directivity = CalculateSignalAtten(Distances, teta, TargetArea, freq, CoordinationFlag);
            directivity = smooth(directivity, 7);
            directivity = directivity'./ max(abs(directivity)); % normailzation so the echo is around 1
            sig_dir = sig_in .* directivity;
            
            % the echo
            sig_echo = multi_echo_gen(sig_dir, n_echoes, time_diff, echo_lvl, snrdB, fs, techo, tend);
            
            % The response
            [FilterBank_echo] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_echo, fs, false, FilterBank_Ref);
            display(['Teta: ' , num2str(teta), ' times: ', num2str(FilterBank_echo.det_ipks/fs,4) ])
            mat_full = reshape([FilterBank_echo.Channels(active_ind).shifted_response], [], n_active);          
            mat_echo = mat_full( ix_mat_echo, : );
            FR_echo = mat_echo(round(nresponse+1), : );
            % take the relevant echo (from detection)
            det_time =  round((techo - tstart)*fs)+1;
            sig_echo = sig_echo( det_time:(det_time + numel(sig_target)) );

            [test_results, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = analyze_echo2ref(sig_target, sig_echo, mat_target, FR_target, mat_echo, FR_echo, ...
                false, pks_target, time_pks_target, f, fs );
            
            test_results.name = string(test_params.name);
            test_results.changed_param = rad2deg(teta); % degrees
            pks_out.name = string(test_params.name);
            pks_out.changed_param = rad2deg(teta); 
            
            %%%% plot %%%%%
            if plot_flags.plot_norm_flag
                my_plot_norm_results(test_results,  sig_target, sig_echo, pks_out, pks_target, time_pks_out, time_pks_target, mat_target, mat_echo, FR_target, FR_echo,   ...
                    mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                    mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2, f, fs);
            end % plot

            %%% update summary table with carreun results
            sum_tbl = [sum_tbl; struct2table(test_results)];
            sum_pks = [sum_pks; pks_out];
            %%% update all Frequecy Reponses
            FR_all_tests(irow,:) = FR_echo;
            irow = irow+1;

        end % for DOA
        
        sum_pks = struct2table(sum_pks);
%         %%%% plot %%%%%
%         if plot_flags.plot_all_reponses
%             my_plot_all_freponses(test_params, FR_target, FR_all_tests, f, fs);
%         end % plot

        %% end DOA Test 
        
    case 'random_noise'
        %%
        detection_num = test_params.tests_num;
        FR_all_tests = nan(detection_num, numel(f));
        irow = 1;
        
        % the echo parameters
        n_echoes = 2;
        echo_lvl = 1;
        noiselvl = 10^(test_params.noise_lvl/10); % changed each test
        techo   = 10e-3;
        tend     = 2e-3;
        % echo- noise pararms
        nsig   = numel(sig_target);
        ndiff  = round(time_diff * fs);
        nend   = round(tend * fs); % nu=
        nstart = round(techo * fs);

        ix_mat_echo = round((nstart-nresponse):(nstart+nresponse));  % the index of the echo respons
        
        
        for k = 1:detection_num

            % the echo - noise 
            sig_echo = noiselvl*randn(1, nstart + nsig + sum(ndiff) + nend);
%             sig_echo = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, snrdB, fs, techo, tend);

            % The response
            [FilterBank_echo] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_echo, fs, false, FilterBank_Ref);
            display(['k: ' , num2str(k), ' times: ', num2str(FilterBank_echo.det_ipks/fs,4) ])
            mat_full = reshape([FilterBank_echo.Channels(active_ind).shifted_response], [], n_active);          
            mat_echo = mat_full( ix_mat_echo, : );
            FR_echo = mat_echo(round(nresponse+1), : );
            % take the relevant echo (from detection)
            det_time =  round((techo - tstart)*fs)+1;
            sig_echo = sig_echo( det_time:(det_time + numel(sig_target)) );

            [test_results, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = analyze_echo2ref(sig_target, sig_echo, mat_target, FR_target, mat_echo, FR_echo, ...
                false, pks_target, time_pks_target, f, fs );
            
            test_results.name = string(test_params.name);
            test_results.changed_param = k;
            pks_out.name = string(test_params.name);
            pks_out.changed_param = k;
            
            %%%% plot %%%%%
            if plot_flags.plot_norm_flag
                my_plot_norm_results(test_results,  sig_target, sig_echo, pks_out, pks_target, time_pks_out, time_pks_target, mat_target, mat_echo, FR_target, FR_echo,   ...
                    mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                    mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2, f, fs);
            end % plot

            %%% update summary table with carreun results
            sum_tbl = [sum_tbl; struct2table(test_results)];
            sum_pks = [sum_pks; pks_out];
            %%% update all Frequecy Reponses
            FR_all_tests(irow,:) = FR_echo;
            irow = irow+1;

        end % for SNR
        
        sum_pks = struct2table(sum_pks);

        %% end random_noise
        
    case 'ChirpBw'
        detection_num = numel(test_params.chirpBW);
        FR_all_tests = nan(detection_num, numel(f));
        irow = 1;
        
        % the echo parameters
        n_echoes = target_params.n_echoes;
        echo_lvl = target_params.echo_lvl;
        snrdBAll    = test_params.SNR;
        techo    = 10e-3;
        tend     = 2e-3;
        
        nstartTarget = tstart*fs;
        nstartEcho = round(techo * fs);
        ix_mat_echo = round((nstartEcho-nresponse):(nstartEcho+nresponse));  % the index of the echo respons
        
        for chirpBW = test_params.chirpBW
            % The Chirp
            fchirp_high = fchirp_low + chirpBW;
            sig_in = PulsePower * chirp(tin, fchirp_high, max(tin), fchirp_low, 'logarithmic'); % chirp(tin, f_high, max(tin), f_low, 'logarithmic')

            %%%%%% THE REFERENCE, and TARGET RESPONSES %%%%% calculated for ecah PULSE
            % reference
            [FilterBank_Ref] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_in, fs,1);
            
            % the target
            sig_target = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, 30, fs, tstart, tend);

            [FilterBank_target] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_target, fs, 0, FilterBank_Ref);
            active_ind = FilterBank_target.Active_Channels;
            n_active = sum(active_ind);
            
            mat_target_full =  reshape([FilterBank_target.Channels(active_ind).shifted_response], [], n_active);

            ix_mat_target = round((nstartTarget-nresponse):(nstartTarget+nresponse));
            mat_target = mat_target_full( ix_mat_target, : );
            ix_start_target = round(nresponse+1);
            FR_target = mat_target(ix_start_target, : );
            f = fc(active_ind);

            % for reference
            [target_results, pks_target, time_pks_target] = analyze_echo2ref(sig_target, sig_target, mat_target, FR_target, mat_target, FR_target, true, [], [], f, fs);
            target_results.name            = "Reference";
            target_results.changed_param   = nan;
            target_results.delay_error     = 0; 
            target_results.is_detected_flg = 0; 
            pks_target.name = "Reference";
            pks_target.changed_param = nan;


            %%%%%% the echo %%%%%%%%%
            for snrdB = snrdBAll
                sig_echo = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, snrdB, fs, techo, tend);
                % The response
                [FilterBank_echo] = FilterBank_Acous_Output(AllFiltersParams, LPF_filt, sig_echo, fs, false, FilterBank_Ref);
                mat_full = reshape([FilterBank_echo.Channels(active_ind).shifted_response], [], n_active);

                % take the relevant echo (from detection)
                actual_echo_time = round((techo - tstart)*fs)+1;
                switch test_params.detection_mode
                    case 'actual'
                        det_time = actual_echo_time;
                        ndelay = 0;
                        is_detected_flg = 1;
                    case 'FB_detection'
                        if any(FilterBank_echo.det_ipks)
                            [k, fb_det] = dsearchn(FilterBank_echo.det_ipks'/fs, (techo - tstart) );
                            det_time = FilterBank_echo.det_ipks(k) - tstart*fs +1;
                            ndelay = det_time - actual_echo_time;
                            is_detected_flg = 1;
                        else
                            det_time = nan;
                            ndelay = det_time - actual_echo_time;
                            is_detected_flg = 0;
                        end % if any(FilterBank_echo.det_ipks)


                end % switch test_params.detection_mode
                % the error of detection
                delay_error = ndelay / fs;

                display(['BW: ' , num2str(chirpBW), ' times: ', num2str(FilterBank_echo.det_ipks/fs,4) , ' delay_error: ', num2str(delay_error,4)])

                if is_detected_flg
                    % the detected signals (shifted)
                    mat_echo = mat_full( (ix_mat_echo + ndelay), : );
                    FR_echo = mat_echo(round(nresponse+1), : );
                    sig_echo = sig_echo( det_time:(det_time + numel(sig_target)-1) );
                    % the result
                    [test_results, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                        mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = analyze_echo2ref(sig_target, sig_echo, mat_target, FR_target, mat_echo, FR_echo, ...
                        false, pks_target, time_pks_target, f, fs );

                else % (not detected) if is_detected_flg
                    % fill nan values
                    sig_echo = sig_echo( actual_echo_time : (actual_echo_time + numel(sig_target)) ); % just for plot
                    [mat_echo, FR_echo, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                        mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2] = deal(nan);
                    test_results = fill_struct_fields(target_results, nan, {'is_detected_flg' , 'name' , 'changed_param', 'delay_error'});
                    pks_out = fill_struct_fields(pks_target, nan, { 'name' , 'changed_param',});
                    time_pks_out = fill_struct_fields(time_pks_target, nan, { 'name' , 'changed_param',});

                end % if is_detected_flg

                test_results.is_detected_flg = is_detected_flg;
                test_results.delay_error     = delay_error;
                test_results.name            = string(test_params.name);
                test_results.changed_param   = chirpBW;
                test_results.snrdB            = snrdB;
                pks_out.name                 = string(test_params.name);
                pks_out.changed_param        = chirpBW;

                %%%% plot %%%%%
                if plot_flags.plot_norm_flag & is_detected_flg
                    my_plot_norm_results(test_results,  sig_target, sig_echo, pks_out, pks_target, time_pks_out, time_pks_target, mat_target, mat_echo, FR_target, FR_echo,   ...
                        mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
                        mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2, f, fs, snrdB);
                end % plot

                %%% update summary table with carreun results
                sum_tbl = [sum_tbl; struct2table(test_results)];
                sum_pks = [sum_pks; pks_out];
                %%% update all Frequecy Reponses
                %             FR_all_tests(irow,:) = FR_echo;
                irow = irow+1;

            end %            for snrdb = snrdBAll

        end % for chirpBW

        sum_pks = struct2table(sum_pks);
%         %%%% plot %%%%%
%         if plot_flags.plot_all_reponses
%             my_plot_all_freponses(test_params, FR_target, FR_all_tests, f, fs);
%         end % plot

        %% end chirp Testt
       

end % switch


%% summary
%%%% Add Chirp Params to sum_table
tt = struct2table(Chirp_params);
tt = repmat(tt, size(sum_tbl,1),1);
sum_tbl = [sum_tbl,tt];

%%%% Summary plot %%%%%
if plot_flags.plot_all_reponses
    norm_flag = 0;
    try
        fig_FR_all = my_plot_all_freponses(test_params, FR_target, FR_all_tests, f, fs);
        fig_pks_analyze = my_plot_pks_tpghs(test_params, sum_pks, sum_tbl );
        fig_all_params = my_params_sumplot(test_params, sum_tbl, norm_flag);
    catch
        warning('plot error')
    end
end % plot

%% save_data
save_flag = 0;
if save_flag
    file_name_prefix = strcat('\', test_params.name, '_') ;
    save_path = 'D:\Dropbox\University\מחקר\experiments\Clutter Influence Simiulation\target_classifier_test\Results_auto';
    
    switch test_params.name
        case 'Interference'
            file_name_prefix =[file_name_prefix, test_params.par2tests, '_'];
        otherwise   %case {'SNR', 'OnePoint_target'}
            add_chirp_params = ['_dur_', num2str(Chirp_params.call_dur*1000), '_fhigh_', num2str(Chirp_params.fchirp_high ), '_flow_', num2str(Chirp_params.fchirp_low)];
            file_name_prefix = [file_name_prefix, add_chirp_params, '_', test_params.detection_mode, '_', 'target_diff_', num2str(target_params.time_diff*1e6)];
    end %  switch test_params.name

    %%%% figures
    
     savefig (fig_FR_all, strcat(save_path, file_name_prefix, 'fig_FR_all') );
    savefig (fig_pks_analyze, strcat(save_path, file_name_prefix, 'fig_pks_alanyze') );
    savefig (fig_all_params, strcat(save_path, file_name_prefix, 'fig_all_params') );

    %%%% tables
    save( strcat(save_path, file_name_prefix, 'tbl_sum_coeff'), 'sum_tbl')
    save( strcat(save_path, file_name_prefix, 'tbl_all_pks'), 'sum_pks')
end % if save_flag
