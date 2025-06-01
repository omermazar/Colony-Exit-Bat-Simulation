function sig_echo = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, snrdB, fs, tstart, tend)

% sig_echo = multi_echo_gen(sig_in, n_echoes, time_diff, echo_lvl, snrdB, fs, tstart, tend)
% returns the ecoh from n_echoes points with time_diff (sec) between them.
% inputs: 
% echo_lvl - the level of the ecgho (linear)
% snrdB - the Singa to Nosise Level (dB)
% fs -sample frequncy
% tstart - the start time of the echo (s)
% tend - quiet time after the end of the echo (s)


nsig   = numel(sig_in);
ndiff  = round(time_diff * fs);
nstart = round(tstart * fs);
nend   = round(tend * fs); % nu=

if n_echoes > 1
    if numel(time_diff)  == 1
        ndiff = [nstart, ndiff * ones(1,n_echoes-1)];
    elseif numel(time_diff) == n_echoes -1
        ndiff = [nstart, ndiff];
    else
        error('Error: The number of delays should be one or n_echoes-1 ')
    end % if time__diff  == 1

    if numel(echo_lvl) == 1
        echo_lvl = echo_lvl * ones(1,n_echoes);
    elseif numel(echo_lvl) ~= n_echoes
        error('Error: The number of echo_lvl should be one or n_echoes ')
    end % if numel(echo_lvl) == 1
else % one point echo
    ndiff = nstart;
end % n_echoes > 2

noiselvl = max(echo_lvl) / 10.^(snrdB/20);
sig_echo = noiselvl*randn(1, nsig + sum(ndiff) + nend + 1);

ix1 = 0;
for k = 1:(n_echoes)
    ix1 = max(2,ix1 + ndiff(k));
    ix = ix1-1 + (1:nsig);
    sig_echo(ix) = sig_echo(ix) + echo_lvl(k) * sig_in;      
end % for k

