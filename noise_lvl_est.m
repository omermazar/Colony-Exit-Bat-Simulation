function [nlvl, nstd] = noise_lvl_est(x,p) 

% estimate thק noise level of a signal aacoording the the lvl in quite periods
% sig_rms = rms(x);
% sig_std = std(x);
% n_ix = abs(x) < sig_rms + sig_std; % 2*sig_std;
% 
% Sep2022: Down samol the signalna for qicker run
if nargin < 2
    p = 70; % 80; % percentile
end
N = numel(x);
% th = prctile(x(1:5:N), p); % With downsamplind
th = prctile(x, p);
n_ix = abs(x) < th;

% Percentile_factor
p_fact = 1.5;
nlvl = p_fact*rms(x(n_ix));
nstd = p_fact*std(x(n_ix));



plot_flag= 0;
if plot_flag
    figure; hold on
    n = 1:numel(x);
    plot(x,'b')
    plot(n(n_ix), x(n_ix), '.r')
    plot(xlim, nlvl*[1 1], '--k')
    plot(xlim, (nlvl+2*nstd)*[1 1], ':k')
end % if plot_flag


