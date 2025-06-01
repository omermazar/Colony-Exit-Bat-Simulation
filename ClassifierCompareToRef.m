function [test_results, pks_out, time_pks_out]  =  ClassifierCompareToRef( FR_target, FR_echo, methods, pks_target, ...
    mat_echo, mat_target, time_pks_target, sig_echo)

%% Compare each peak reasponse to the refence rsponse of the desired object (target)
%   
%   test_results = ( ClassifierCompareToRef( FR_target, FR_echo, methods, pks_target) 
%               returns a struct with the results of the tests defindes by
%               methods.
%
%   [test_results, pks_out, time_pks_out]  =  ClassifierCompareToRef(___)
%       returns also structs with the parameterd of the peaks in the frequency response and time, respectivelty
%
% [test_results, pks_out, time_pks_out]  =  ClassifierCompareToRef( FR_target, FR_echo, methods, pks_target, ...
%       mat_echo, mat_target, time_pks_target, sig_echo)
%
%  Inputs: 
%       FR_target: the frequency response of the target
%       FR_echo:  the frequency response of the target                  
%       methods: methods a struct of the desired methods to aply on classifier, wtih the following default parametrs :
%           methods.target_flg        = 0;  % for refence parmeters define as 1
%           methods.normalzation      = 1;  % noramlize the signals according
%                                       % to their energy 
%           methods.standarzization   = 0;  % standardization of the singnals
%                                           % (mean =1; std =0), (Doesnt work well)
%                                           % Should be standarzization or  normalzation                           
%           methods.correaltionTest   = 1;  % test correlation coefficent between the signals% 
%           methods.pksTrghsTest      = 1;  % test correlation coefficent between the signals% 
%           methods.mseTest           = 1;  % test correlation coefficent between the signals% 
%           methods.cochleogramMatrix = 0;  % run calssifier on the time-frequency response matrix
%           methods.frequencyResponse = 1;  % run calssifier on the frequency response vector
%           methods.timeResponse      = 0;  % run calssifier on the time response vector
%       pks_target: a struct with peaks of the refernce teraget
%       mat_echo: the relevant cochlogram of the echo, required only if
%           methods.cochleogramMatrix == 1
%       mat_target: the relevant cochlogram of the echo (for reference)
%       time_pks_target: a sstuct of the detected times od the target
%           signal, required only if  methods.timeResponse == 1
%       sig_echo: the signal itself required only if  methods.timeResponse == 1

% Based on analyze_echo2ref()

%% Input
target_flg = methods.target_flg;
% if nargin < 4
%     pks_target = [];
% end
% if nargin < 5
%     mat_echo = nan;
% end
% if nargin < 6
%    mat_target = nan;
% end % if nargin <6 
% 
% if nargin < 7
%    time_pks_target = [];
% end % if nargin <6 
% 
% if nargin < 8
%    sig_echo = nan;
% end % if nargin <6 


%% normailzation
% Standardiazion - 'Kills the signals'
if methods.standarzization
    mat_target = normalize(mat_target);
    FR_target  = normalize(FR_target);
    mat_echo   = normalize(mat_echo);
    FR_echo    = normalize(FR_echo);
end % if methods.standarzization

% Noramazilation equal energy
if methods.normalzation
    FR_target  = FR_target ./ norm(FR_target);
    FR_echo    = FR_echo ./ norm(FR_echo);
end % if methods.normalzation

%% Peaks and troughs Test
% Remark: Thresholds were testet for 4ms 75-70 ms

if methods.pksTrghsTest
    % equal energy
%     if methods.normalzation
%         FR_target  = FR_target ./ norm(FR_target);
%         FR_echo    = FR_echo ./ norm(FR_echo);
%     end % if methods.normalzation

    if target_flg
        pks_target = [];
        %%% refernce
        pr_TH = 0.02; % 0.05;
        pk_Th = 0.02; % 0.05;
        min_dist = 3 ; %3e3; % Hz

    else % ~target_flg
        pr_TH =  0.02; % 0.05; %0.02;
        pk_Th =  0.02; % 0.05; %0.02;
        min_dist = 3 ;% 3e3; % Hz

    end  %

    %%% Freq- peaks
    [pks_out, pks_coeff] = pks_trghs_coeff(FR_echo, pk_Th, pr_TH, min_dist,  target_flg, pks_target, false);

    %%%% output
    test_results.pks_trghs_rms = pks_coeff.pks_trghs_rms;
    test_results.pks_rms       = pks_coeff.pks_rms;
    test_results.trghs_rms     = pks_coeff.trghs_rms;
end % pksTrghsTest


%% correlation test 
if methods.correaltionTest && ~target_flg

    if methods.frequencyResponse
        % Correlation Max Value
%         corr = xcorr(FR_target, FR_echo);
%         test_results.correlationMaxFR = max(corr);

        % Correlation Coefficient
        FR_cor_coeff = corrcoef(FR_target, FR_echo);
        test_results.corrCoeffFR = FR_cor_coeff(2,1);
    end % if frequencyResponse

    if methods.cochleogramMatrix
        % Noramziaion
        if methods.normalzation
            mat_target = mat_target ./ norm(mat_target);
            mat_echo   = mat_echo ./ norm(mat_echo);
        end % if methods.normalzatio
        
         % Correlation Max Value
%         corr = xcorr2(mat_target, mat_echo);
%         test_results.correlationMaxMatrix = max(corr, [], 'all')';
%         [corr_pks, row, column] = my_findpeakds2d(corr);

        % Correlation Coefficient
        mat_corr_coeff = corrcoef(mat_target, mat_echo);
        test_results.corrCoeffMatrix = mat_corr_coeff(2,1);
        
        % Correlation Shift in Time
        % three most dominant peaks
%         n_pks = min(3,numel(corr_pks));
%         if n_pks > 0
%             [p,ii] = maxk(corr_pks,n_pks);
%             corr_ix = row(ii);
%             corr_shift = mean(corr_ix) - (size(mat_echo,1) +1); % the middle of the correlation rows is zero
%         else
%             corr_shift = nan;
%         end
%         %         [m, mi] = max(corr2);
%         test_results.correlationShiftMatrix = round(corr_shift);
    end % methods.frequencyResponse


end % if  methods.correaltionTest

%% Normalized Distance between echoes (MSE), Similarity Test
if methods.mseTest && ~target_flg
    if methods.frequencyResponse
        test_results.mseFR = mse(FR_target, FR_echo);
    end % if methods.frequencyResponse

    if methods.cochleogramMatrix
        test_results.mseMatrix = mse(mat_target, mat_echo);
        test_results.ssimMarix = ssim(mat_target, mat_echo);
    end % if methods.cochleogramMatrix
end % if methods.mseTest


%% time analysis
if methods.timeResponse
    if methods.normalzation
        sig_echo   = sig_echo ./ norm(sig_echo);
    end % if methods.normalzatio
    % findpeaks params
    t_pr_TH = 0.025;
    t_pk_Th = 0.025;
    t_min_dist = 10 ; %50 micros; %

    env1 = smooth(envelope(sig_echo),5)';
    %%%% time peaks
    [pks_detected, pks_coeff] = pks_trghs_coeff(env1, t_pk_Th, t_pr_TH, t_min_dist,  target_flg, time_pks_target , true);

    %%%% output
    test_results.time_pks_trghs_rms = pks_coeff.pks_trghs_rms;
    test_results.time_pks_rms       = pks_coeff.pks_rms;
    test_results.time_trghs_rms     = pks_coeff.trghs_rms;

    time_pks_out = pks_detected;
end % methods.timeResponse 

