function [detected_echoes] = ...
    FilterBank_echoes_decision( ...
    fc_detect_struct, fb_ratio_th, delay_tollerance, ...
    CurrEchosFromPreyStruct, DetectedPreysVec, rel_InterferenceFullStruct, FilterBank, fs_ts, det_method)
% function [] = FilterBank_echoes_decision()
% this function integrates the detections of each fc and determines the detected ehhoes of the current pulse
% the function also tries to link the detected echoes to targets

%% init all output
detected_delays = [];
all_errs = [];
linked_delays = [];
estimated_targets = [];
est_errs = [];
unlinked_delays = [];
ul_jammer = [];
detected_counts = [];
detections_idx  = [];
fb_detect_th = [];
detected_width = [];
all_rise_times  = [];
det_time_ind  = [];
powers_db  = [];
id_small_errs=[];
id_linked = [];

%% General inputs
all_delays = [fc_detect_struct.detection_delays];
all_rise_times = [fc_detect_struct.rise_time_det];
all_delays = all_delays(all_delays>0);
all_rise_times = all_rise_times(all_delays>0);
all_powers = [fc_detect_struct.rx_power_db];

it_max_num = 10;
% % [predetected_delays, detected_counts, detected_width, cluster_ind] = ...
% %     cluster_by_dist(all_delays, delay_tollerance );
fcs_active_num = numel(fc_detect_struct);
fb_detect_th = floor(fb_ratio_th *fcs_active_num);
%%

switch det_method
    
    case 'Integrated'
        %% Integrated detection
        detected_delays =  all_delays;
        powers_db = all_powers;
        det_time_ind = true(size(all_rise_times));
        
    case 'Fcs'
        %% FC's detection
        %% Cluster the delays by maximum delay-tolleance
        
        if ~isempty(all_delays)
            
            [clustered_delays_strct] = ...
                cluster_by_dist(all_delays, delay_tollerance, it_max_num );
            
            detected_counts = [clustered_delays_strct.counts];
            detected_delays_temp = [clustered_delays_strct.data];
            detected_width = [clustered_delays_strct.width];
            
            % the rise times of ech detection of the echoes
            det_time_ind = vertcat(clustered_delays_strct.data_ind)';% the index of the raw data input to clusters
            
            detections_idx = detected_counts >= fb_detect_th;
            detected_delays = detected_delays_temp(detections_idx);
            
            %% calculate powers of detection
            powers_db = zeros(1,numel(clustered_delays_strct));
            for k = 1:numel(clustered_delays_strct)
                if detections_idx(k)
                    idk = clustered_delays_strct(k).all_ind;  % idk = clustered_delays_strct(k).data_ind;
                    powers_db(k)= mean(all_powers(idk));
                end % if
            end % for

        end  %  if ~isempty(all_delays)
end % switch det_method

%% find the detected target according to FilterBank_Targets_track

if ~isempty(detected_delays)
    %% compare detection to 'real targets  ' - referene, errors etc.
    echos_start_times = CurrEchosFromPreyStruct.EchosTimes*fs_ts; % compare to tx time
    [detected_prey_idx, d] = dsearchn(CurrEchosFromPreyStruct.TargetIndex', DetectedPreysVec');
    detected_prey_idx = detected_prey_idx(d==0);
    
    detected_echos_st_ref = round(echos_start_times(detected_prey_idx) );
    
    if ~isempty(detected_echos_st_ref) % XXX
        [vi, ~] = dsearchn(detected_echos_st_ref', detected_delays' );
        all_errs = detected_delays- detected_echos_st_ref(vi');
        
        id_small_errs = abs(all_errs) <= delay_tollerance;
        f_small = find(id_small_errs);
        %     estimated_targets = zeros(size(detected_delays));
        estimated_targets = DetectedPreysVec(vi(id_small_errs));
        [estimated_targets, ia, ~] = unique(estimated_targets,'stable');
        
        est_errs = all_errs(id_small_errs);
        
        est_errs = est_errs(ia);
        id_linked = zeros(size(id_small_errs));
        id_linked(f_small(ia)) = id_small_errs(f_small(ia));
        id_linked = logical(id_linked);
        
        linked_delays = detected_delays(id_linked);
        [unlinked_delays, idx_unlinked] = setdiff(detected_delays, linked_delays);
        
       
    else % if ~isempty(detected_echos_st_ref)
        %             estimated_targets = [];
        
        %             det_time_ind =[];
    end % if ~isempty(detected_echos_st_ref)
else % if ~iseampty(detected_delays)
    %         estimated_targets = [];
    
end %if ~iseampty(detected_delays)

% % %% calculate powers of detection
% % powers_db = zeros(1,numel(clustered_delays_strct));
% % for k = 1:numel(clustered_delays_strct)
% %     if detections_idx(k)
% %         idk = clustered_delays_strct(k).all_ind;  % idk = clustered_delays_strct(k).data_ind;
% %         powers_db(k)= mean(all_powers(idk));
% %     end % if
% % end % for

%% try to find jamming bat causing the unlinked
% detected_fcs_idx= horzcat(fc_detect_struct.active_fc_idx);
% mean_delay_ir = mean(FilterBank.filter_delay_times(detected_fcs_idx));
% std_delays_ir =  std(FilterBank.filter_delay_times(detected_fcs_idx));
% if ~isempty(unlinked_delays)
%     rel_masking_times = vertcat(rel_InterferenceFullStruct.Times_precise)*fs_ts;
%     rel_jammers = vertcat(rel_InterferenceFullStruct.BATTxNum); % the bat-tx in that relevant times
%     
%     ul_jammer = zeros(size(unlinked_delays));
%     if ~isempty(rel_masking_times)
%         
%         for kk = numel(idx_unlinked)
%             times_incluster = all_rise_times(clustered_delays_strct(idx_unlinked(kk)).all_ind)';
%             [closest_jammer_idx, d] = dsearchn(rel_masking_times, (times_incluster-mean_delay_ir) );
%             th_idx = closest_jammer_idx(d < 2*(std_delays_ir+delay_tollerance) ); % 2*delay_tollerance
%             if ~isempty(th_idx)
%                 jammers = rel_jammers(th_idx);
%                 if numel(jammers) > 1
%                     u_jammers = unique(jammers);
%                     jammer_counts = nonzeros(accumarray(jammers,jammers',[],@numel));
%                     [~, imax] = max(jammer_counts);
%                     ul_jammer(kk) = u_jammers(imax);
%                 else % if numel(jammers) > 1
%                     ul_jammer(kk) = jammers;
%                 end % if numel(jammers) > 1
%             end
%         end % for kk
%     end % if isempty(rel_masking_times)
% end % if ~isempty(unlinked_delays





%% echoes detected- output
detected_echoes = struct( ...
    'all_detected_delays', round(detected_delays(:))', ...
    'delays_samples_errors', all_errs(:)', ...
    'linked_delays', round(linked_delays(:))', ...
    'linked_idx', id_linked, ...        %   ,id_small_errs
    'estimated_targets', estimated_targets(:)', ...
    'estimated_errors', est_errs(:)', ...
    'unlinked_delays', round(unlinked_delays(:))', ...
    'unlinked_est_jammer', ul_jammer(:)', ...
    'cluster_counts', detected_counts(detections_idx), ...
    'counts_th', fb_detect_th, ...
    'cluster_width', detected_width(detections_idx), ...
    'all_detected_rise_times', round(all_rise_times(det_time_ind)), ...
    'powers_db', nonzeros(powers_db)', ...
    'link_th', delay_tollerance);



%% analyze plot
plot_flag = 0;
if plot_flag
    plot_echoes_histogram
end