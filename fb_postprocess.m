function [post_process_mat] = fb_postprocess(curr_FilterBank_Pulses_strct, rel_rx_mat, FilterBank, pp_method, pp_gain)

% inputs: 
%       curr_FilterBank_Pulses_strct - BAT.FilterBank_Pulses_strct
%         rel_rx_mat- the rx mat in each fc in the relevant times
%         FilterBank - the sturct of FB for this run
%         pp_method = 4- 'lateral inhibition', 0 - 'integrated mode only'  

%% Test corssTalk

% curr_FilterBank_Pulses_strct = BAT(BatNum).FilterBank_Pulses_strct(CurrentPulseNum);
u_fcs_idx = unique(curr_FilterBank_Pulses_strct.TxActiveFcs_idx, 'stable');
ct_mat = FilterBank.cross_matrix(u_fcs_idx, u_fcs_idx);
n_rows = size(ct_mat,1);
% ct_1 = inv(ct_mat);


if pp_method ==0
    post_process_mat = [];

elseif pp_method == 1 % 
    ct_1 = -ct_mat;
    diag_ct = sub2ind(size(ct_1), 1:n_rows, 1:n_rows);
    ct_1(diag_ct) = n_rows/2*ones(1,n_rows); % n_rows
    % a_norm = diag(1./sum(ct_1')); 
    k1= 2;
    post_process_mat = k1*ct_1*rel_rx_mat(u_fcs_idx,:);
    
elseif pp_method ==2
    ct_1 = eye(n_rows) - 1/n_rows*ones(n_rows);
    k1= 1;
    post_process_mat = k1*ct_1*rel_rx_mat(u_fcs_idx,:);
    
elseif pp_method ==3
    temp_mat= ones(n_rows);
    % the neaest chennle id much reduced
    diag1_ind= [(n_rows+1):(n_rows+1):numel(temp_mat)];
    diag2_ind= [2:(n_rows+1):numel(temp_mat)];
%     diag3_ind= [(n_rows+2):(n_rows+1):numel(temp_mat)];
%     diag4_ind= [3:(n_rows+1):numel(temp_mat)];
    k1= 3;
%     k2= 2;
    temp_mat([diag1_ind, diag2_ind]) = k1;
%     temp_mat([diag3_ind, diag4_ind]) = k2;
    ct_1 = k1*eye(n_rows) - 1/n_rows*temp_mat;
  
    post_process_mat = k1*ct_1*rel_rx_mat(u_fcs_idx,:);

% elseif pp_method == 5 % delay time
% % %      n_rows = size(ct_mat,1);
% % %      % each time the freq changes
% % %      tx_fcs_start_times = [0, find(diff(curr_FilterBank_Pulses_strct.TxActiveFcs_idx))];
% % %      ct_1 = 2*eye(n_rows)-1*ones(n_rows);

elseif pp_method == 4 % matrices
    ct_mat = round(ct_mat,2);
    ct_1 = ct_mat';
    A = (ct_1*ct_mat)\ct_1;
%     dA = diag(A);
%     A = 0.5*A; 
%     A(1:n_rows+1:end) = abs(dA);
    post_process_mat = pp_gain*A*rel_rx_mat(u_fcs_idx,:);
%     post_process_mat = (ct_1*ct_mat)\ct_1* rel_rx_mat(u_fcs_idx,:);
%     post_process_mat = inv(ct_1*ct_mat) * ct_1* rel_rx_mat(u_fcs_idx,:);
end % if pp_method

%% Integratet channel results

% in_mat = rel_rx_mat(u_fcs_idx,:);  % post_process_mat; % rel_rx_mat(u_fcs_idx,:)
% % delay_mat = zeros(size(rel_rx_mat(u_fcs_idx,:)));
% delay_mat = zeros(size(in_mat));
% % sm_mat = delay_mat; % XXXXX
% % sorted_fcx = sort(u_fcs_idx,'descend'); 
% tx_delays = [0, find(diff(curr_FilterBank_Pulses_strct.TxActiveFcs_idx))'];
% % in_mat = ref_mat; % rel_rx_mat(u_fcs_idx,:);
% 
% sm_coeff = mean(tx_delays/5); % 10 - too narrow , 2 - too wide
% 
% % % figure
% for k = 1:numel(u_fcs_idx)
%     d = tx_delays(k);
%     zero_pad = zeros(1,d);
%     delay_mat(k,:) = smooth([ in_mat(k,d+1:end), zero_pad], sm_coeff);
% %     sm_mat(k,:) = smooth(delay_mat(k,:), sm_coeff);
% 
% %     plot(in_mat(k,:))
% % %     figure
% % %     hold on
% % %     plot(delay_mat(k,:))
% % %     plot(sm_mat(k,:));
% % %     title(num2str(u_fcs_idx(k)))
% end % for k
% int_channel = sum(delay_mat); % *1./numel(u_fcs_idx);
% % sm_ch = sum(sm_mat);  % XXXX
% % int_channel = sm_ch;  % XXXX


% 