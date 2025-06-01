function [int_ch ] = fb_integrate_filters(rx_struct, rel_rx_mat)

% summation of all the rx-filters with delay lines according to the tx
% pulse
% inputs:   rx_struct = the strauct of fb rx powers and times
%           rel_rx_mat = the matrix after convolution - only in the tx
%           relevant fc's

% picking the active frequnciues in the vurrent pulse
u_fcs_idx = unique(rx_struct.TxActiveFcs_idx, 'stable');
in_mat = rel_rx_mat;
% in_mat = rel_rx_mat(u_fcs_idx,:);  % post_process_mat; % rel_rx_mat(u_fcs_idx,:)

% calc the delay of each tx frequency-band (fc)
tx_delays = [0, find(diff(rx_struct.TxActiveFcs_idx))'];

% smoothing parameter (lpf of each fc -channel)
sm_coeff = mean(tx_delays/5); % 10 - too narrow , 2 - too wide

%init output
delay_mat = zeros(size(in_mat));

% delay and zero-pad each channel
% % figure
for k = 1:numel(u_fcs_idx)
    d = tx_delays(k);
    zero_pad = zeros(1,d);
    delay_mat(k,:) = smooth([ in_mat(k,d+1:end), zero_pad], sm_coeff);
%     sm_mat(k,:) = smooth(delay_mat(k,:), sm_coeff);

%     plot(in_mat(k,:))
% %     figure
% %     hold on
% %     plot(delay_mat(k,:))
% %     plot(sm_mat(k,:));
% %     title(num2str(u_fcs_idx(k)))
end % for k
int_ch = sum(delay_mat); % *1./numel(u_fcs_idx);


