function [mat_in, u_fft_idx ] = calc_jam_raw_mat(mat_in, orig_time_vec, freq_in, power_in, start_us_time, fft_freqs, fs_ts)

% mat_out= mat_in;
% mat_out= zeros(size(mat_in));

% temp_mat = zeros(size(mat_in));
max_samples = size(mat_in,2);


ref_time_vec_us = linspace(orig_time_vec(1)*fs_ts,  (orig_time_vec(end)+1)*fs_ts, ...
        max(2,numel(freq_in)) );
time_vec_us = ref_time_vec_us(1) : ref_time_vec_us(end)-1;
time_vec_idx = round(time_vec_us - start_us_time); % rows


% protect end of signal
vec_terminate = time_vec_idx <= max_samples ;

if ~isempty(vec_terminate) && max_samples > 1
    time_vec_us = time_vec_us(vec_terminate);
    time_vec_idx = time_vec_idx(vec_terminate);
    
    freqs_us = interp1(ref_time_vec_us, freq_in, time_vec_us)*1000; % ?Hz
    power_us = interp1(ref_time_vec_us, power_in, time_vec_us);
    max_value = max(fft_freqs);
    if max(freqs_us) > max_value
        freqs_us(freqs_us>= max_value) = max_value;
    end
%     mat_ind = zeros(size(power_us));
    
    v_fft = 1:numel(fft_freqs);
   
    fft_idx = round(interp1(fft_freqs, v_fft, freqs_us)); % columes
    u_fft_idx = unique(fft_idx);
    
    % update ooutput matrix
    mat_ind = sub2ind(size(mat_in), fft_idx, time_vec_idx  );

%     mat_out(mat_ind) = mat_out(mat_ind) + power_us;
%     vec_out = mat_out(mat_ind) + power_us; 
%     vec_out = mat_in(mat_ind) + power_us;
    mat_in(mat_ind) = mat_in(mat_ind) + power_us; % vec_out;
%     mat_out(mat_ind) = power_us;

% %     mat_out= mat_out + temp_mat;
        
%     mat_out(mat_ind) = mat_out(mat_ind) + power_us;

else % ~isempty(vec_terminate)
    u_fft_idx = [];
%     mat_out = mat_in;
end % if ~isempty(vec_terminate)