function [test_results, pks_out, time_pks_out, mat_target_n1, mat_echo_n1, FR_target_n1, FR_echo_n1, ...
    mat_target_n2, mat_echo_n2, FR_target_n2, FR_echo_n2]  = analyze_echo2ref(sig_target, sig_echo, mat_target, FR_target, mat_echo, FR_echo, ...
    target_flg, pks_target, time_pks_target, f, fs)

%% normailzation
% Standardiazion - 'Kills the signals'
n_FR = numel(FR_target);
mat_target_n1 = normalize(mat_target);
FR_target_n1 = normalize(FR_target);
mat_echo_n1 = normalize(mat_echo);
FR_echo_n1 = normalize(FR_echo);

% equal energy
mat_target_n2 = mat_target ./ norm(mat_target);
FR_target_n2 = FR_target ./ norm(FR_target);
mat_echo_n2 = mat_echo ./ norm(mat_echo);
FR_echo_n2 = FR_echo ./ norm(FR_echo);

sig_target = sig_target ./ norm(sig_target);
sig_echo = sig_echo ./ norm(sig_echo);

%% correlation test 
corr = xcorr(FR_target_n1, FR_echo_n1,'biased');
test_results.FR_corr_n1 = max(corr);

corr = xcorr(FR_target_n2, FR_echo_n2);
test_results.FR_corr_n2 = max(corr);

corr = xcorr2(mat_target_n1, mat_echo_n1);
test_results.mat_corr_n1 = max(corr, [], 'all')';


corr = xcorr2(mat_target_n2, mat_echo_n2);
test_results.mat_corr_n2 = max(corr, [], 'all')';

%%% find the time diff between the signals
% [corr2, l2] = xcorr(sum(mat_target_n2,2), sum(mat_echo_n2,2));
% [m, mi] = max(corr2);
% test_results.corr_shift = l2(mi);
corr2 = xcorr2(mat_target_n2, mat_echo_n2 ); % 2D correlaetiob
[corr_pks, row, column] = my_findpeakds2d(corr2);
% three most dominant peaks
n_pks = min(3,numel(corr_pks));
if n_pks > 0
    [p,ii] = maxk(corr_pks,n_pks);
    corr_ix = row(ii);
    corr_shift = mean(corr_ix) - (size(mat_echo,1) +1); % the middle of the correlation rows is zero
else
    corr_shift = nan;
end
[m, mi] = max(corr2);
test_results.corr_shift = round(corr_shift);

mat_corr_coeff = corrcoef(FR_target, FR_echo);
test_results.mat_corr_coeff = mat_corr_coeff(2,1);

FR_cor_coeff = corrcoef(FR_target, FR_echo);
test_results.FR_cor_coeff = FR_cor_coeff(2,1);

%% Normalized Distance between echoes (MSE)
test_results.FR_mse_n1 = mse(FR_target_n1, FR_echo_n1);
test_results.FR_mse_n2 = mse(FR_target_n2, FR_echo_n2);

test_results.mat_mse_n1 = mse(mat_target_n1, mat_echo_n1);
test_results.mat_mse_n2 = mse(mat_target_n2, mat_echo_n2);

%% Structural similarity 
test_results.mat_ssim_n2 = ssim(mat_target_n2, mat_echo_n2);
test_results.mat_ssim_n1 = ssim(mat_target_n1, mat_echo_n1);

%% Peaks and troughs
%%% refernce
% Thresholds were testet for 4ms 75-70 ms
if target_flg
    pr_TH = 0.02; % 0.05;
    pk_Th = 0.02; % 0.05;
    min_dist = 3 ; %3e3; % Hz

else % ~target_flg
    pr_TH =  0.02; % 0.05; %0.02;
    pk_Th =  0.02; % 0.05; %0.02;
    min_dist = 3 ;% 3e3; % Hz
    
    corr_pks = pks_target.pks;
    ipks = pks_target.ipks;
    prs_pks = pks_target.prs_pks;
    trghs = pks_target.trghs;
    itrghs = pks_target.itrghs;
    prs_trghs = pks_target.prs_trghs;
end  %

% [pks, ipks, ~, prs_pks] = findpeaks(FR_target_n2, 'MinPeakHeight',pk_Th, 'MinPeakProminence', pr_TH, 'MinPeakDistance', min_dist);
% [trghs, itrghs, ~, prs_trghs] = findpeaks(-FR_target_n2, 'MinPeakProminence', pr_TH, 'MinPeakDistance', min_dist);

%%% Freq- peaks
[pks_detected, pks_coeff] = pks_trghs_coeff(FR_echo_n2, pk_Th, pr_TH, min_dist,  target_flg, pks_target, false);

%%%% output
test_results.pks_trghs_rms = pks_coeff.pks_trghs_rms;
test_results.pks_rms       = pks_coeff.pks_rms;
test_results.trghs_rms     = pks_coeff.trghs_rms;


pks_out = pks_detected;

%% time analysis
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

