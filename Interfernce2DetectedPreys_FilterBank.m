 function [FindsMaskingStruct ] = Interfernce2DetectedPreys_FilterBank( CurrPulseNum, ...
    FilterBank_nomask_echoes, FilterBank_all_echoes, ...
    CurrEchosFromPreyStruct, InterferenceFullStruct, BAT, BatNumRx, AllParams)


% function [FindsMaskingStruct ] = Interfernce2DetectedPreys_FilterBank( ...
%     FilterBank_detected_echoes, FilterBank_jamming_echoes, ...
%     CurrEchosFromPreyStruct, InterferenceFullStruct, CurrPulseNum, ...
%     BAT, BatNum, AllParams, FilterBank )

% %Interfernce2DetectedPreys_FilterBank() returns a sturct of the masking to the
% detected prey items including the loclization errors
%  The function is for FilterBank Mode , replacing Interfernce2DetectedPreys() for this mode


%% decide whether false alarmes are masking or inside 'filter noise' 
time_link_th = FilterBank_nomask_echoes.link_th;
power_link_th = 1; % larger diffence is not the same signal db
filterbank_fs_ts =  AllParams.BatSonarParams.FilterBank_Fs * AllParams.SimParams.SampleTime;


NumOfSamples = AllParams.SimParams.SimulationTime / AllParams.SimParams.SampleTime + 1;

det_prey_nomaksing = FilterBank_nomask_echoes.estimated_targets;
masked_prey_idx = false(size(det_prey_nomaksing));

if isempty(FilterBank_nomask_echoes.all_detected_delays) % no detections withoutmask
    idx_jammer = ones(size(FilterBank_all_echoes.all_detected_delays));

else % if isempty(FilterBank_nomask_echoes.all_detected_delays)
    
    tx_time = FilterBank_nomask_echoes.all_detected_rise_times(1) - FilterBank_nomask_echoes.all_detected_delays(1);
    
    [id_closest_nm, dtime] = dsearchn(FilterBank_nomask_echoes.all_detected_delays', ...
        FilterBank_all_echoes.all_detected_delays');
    id_closest_nm = id_closest_nm(:)';
    dtime = dtime(:)';
    
    time_jam_idx = dtime > time_link_th;
    dpower = FilterBank_all_echoes.powers_db - FilterBank_nomask_echoes.powers_db(id_closest_nm);
    power_jam_idx = dpower > power_link_th;
    
    idx_jammer = time_jam_idx(:)' | power_jam_idx(:)';
    
end % if isempty(FilterBank_nomask_echoes.all_detected_delays)


if sum(idx_jammer) > 0 % if thre are jammers
    
    %% decide wether the target-echo is masked
    linked_not_masked_idx =  FilterBank_all_echoes.linked_idx & ~idx_jammer;
    unmasked_linked_delays = FilterBank_all_echoes.all_detected_delays(linked_not_masked_idx);
    [~, um_targets_idx] = intersect(FilterBank_all_echoes.linked_delays, unmasked_linked_delays);
    
    unmasked_prey = FilterBank_all_echoes.estimated_targets(um_targets_idx);
    
    [masked_prey, i_mp] =  setdiff(FilterBank_nomask_echoes.estimated_targets, unmasked_prey, 'stable');
%     masked_prey_idx = false(size(det_prey_nomaksing));
    masked_prey_idx(i_mp) = 1;
    IsAnyPreyMasked = ~isempty(masked_prey);
    
    %% try to link tx-masking bat to jamming
    jammer_times = FilterBank_all_echoes.all_detected_rise_times(idx_jammer);
    jammer_delays = FilterBank_all_echoes.all_detected_delays(idx_jammer);
    jammer_powers = FilterBank_all_echoes.powers_db(idx_jammer);
    
    estimated_tx_masking_bat = zeros(size(jammer_times));
    
   
    ds_time_vec = [round((min(jammer_times)-time_link_th)/filterbank_fs_ts) :...
        round((max(jammer_times)+time_link_th)/ filterbank_fs_ts)];
    rel_Interference_strct =  BAT(BatNumRx).InterferenceFullStruct(ds_time_vec);
    all_int_times = vertcat(rel_Interference_strct.Times_precise) * filterbank_fs_ts;
    all_int_BatNum = vertcat(rel_Interference_strct.BATTxNum);
    
    if ~isempty(all_int_times)
        [k_times, time_diff] = dsearchn(all_int_times, jammer_times');
        k_times = k_times(:)';
        time_diff = time_diff(:)';
        
        ia = time_diff <= (time_link_th+filterbank_fs_ts);
        k_int = k_times(ia);
        estimated_tx_masking_bat(ia) = [all_int_BatNum(k_int)];
    end % if ~isempty(all_int_times)
    
    if IsAnyPreyMasked
        %     [masked_prey_echo_idx, d] = dsearchn(CurrEchosFromPreyStruct.TargetIndex', masked_prey');
        %     masked_prey_echo_idx = masked_prey_echo_idx(d==0);
        masking_delay = FilterBank_nomask_echoes.linked_delays(i_mp);
        
        TotalMaskingTimes = round((tx_time + masking_delay) / filterbank_fs_ts);
         
    else
        TotalMaskingTimes = [];
    end % if IsAnyPreyMasked
else % if sum(idx_jammer > 0)
   
    %% there is no jamming 
    masked_prey = [];
    unmasked_prey = det_prey_nomaksing;
    unmasked_linked_delays = FilterBank_nomask_echoes.linked_delays;
    
    IsAnyPreyMasked = 0;
    TotalMaskingTimes = [];
    jammer_delays = [];
    jammer_powers = [];
    estimated_tx_masking_bat = [];

end % if sum(idx_jammer > 0)
%% SIR calculation
% for each estmiated target find the SIR to closest Interferef if it is
% in the masking time window % less than 2*TH
SI_bwd_time =  AllParams.BatSonarParams.MaskingBckdTime * AllParams.BatSonarParams.FilterBank_Fs*1e-3;
SI_fwd_time =  AllParams.BatSonarParams.MaskingFwdTime * AllParams.BatSonarParams.FilterBank_Fs*1e-3;

SIR_to_detected = FilterBank_nomask_echoes.powers_db(FilterBank_nomask_echoes.linked_idx); % start with rx power
% for each detected find the relevant making signal
if ~isempty(jammer_delays)
    for m = 1:numel(FilterBank_nomask_echoes.linked_delays)
        cur_delay= FilterBank_nomask_echoes.linked_delays(m);
        min_period = cur_delay - SI_fwd_time;
        max_period = cur_delay + SI_bwd_time;
        rel_jam_id = ( jammer_delays>=min_period & jammer_delays<=max_period );
        
        if sum(rel_jam_id)>0
            max_jam = max(jammer_powers(rel_jam_id));
            if max_jam > AllParams.BatSonarParams.NoiseLeveldB
                SIR_to_detected(m) = SIR_to_detected(m) - max_jam;
            end  % if max_jam > AllParams.BatSonarParams.NoiseLeveldB
            
        end %if ~isempty(rel_jam)

    end %for m
    
end % if jammer_delays

% % if ~isempty(jammer_delays)
% %     [id_closest_masker, d2masker] = dsearchn(jammer_delays', FilterBank_nomask_echoes.linked_delays');
% %     id_closest_masker = id_closest_masker(:)';
% %     d2masker = d2masker(:)';
% %     masked_idx = find(d2masker < 2*time_link_th);
% %     
% %      if ~isempty(masked_idx)
% %         SIR_2all = FilterBank_nomask_echoes.powers_db(FilterBank_nomask_echoes.linked_idx) - jammer_powers(id_closest_masker);
% %         SIR_to_detected(masked_idx) = min([SIR_to_detected(masked_idx) ; SIR_2all(masked_idx)]);
% %     end % if ~isempty(masked_idx)
% % end % if jammer_delays




%% localization errors

%%% errors to target detections
% no-masking targets (a;; target withput masking)

nm_prey_powerdb = FilterBank_nomask_echoes.powers_db(FilterBank_nomask_echoes.linked_idx);
nm_prey_times = FilterBank_nomask_echoes.all_detected_rise_times(FilterBank_nomask_echoes.linked_idx);

% masked targets - SIR =0 db
DetecetedPrey2InterferenceRatioDB = SIR_to_detected;
% DetecetedPrey2InterferenceRatioDB(i_mp) = 0;

if AllParams.BatSonarParams.LocalizationErrFlag
    [DFErr, RangeErr, RelativeDirectionErr] = CalculateLocalizationErrors(...
        BAT(BatNumRx), CurrEchosFromPreyStruct, det_prey_nomaksing, ...
        DetecetedPrey2InterferenceRatioDB, nm_prey_powerdb, CurrPulseNum, AllParams);
else % if LocalizationErrFlag
    DFErr =0;
    RangeErr = 0;
    RelativeDirectionErr = 0;
end % if LocalizationErrFlag

%% find masking interferences frequeincies for JAR response

if IsAnyPreyMasked
    interval_edge = AllParams.BatSonarParams.PulseDuration_Search*1e-3/AllParams.SimParams.SampleTime;
    InterMinFreq = zeros(size(TotalMaskingTimes));
    for m=1:numel(TotalMaskingTimes)
        rel_vec_idx = [max((TotalMaskingTimes(m)- ceil(interval_edge/2)),1):...
            min((TotalMaskingTimes(m) + interval_edge),NumOfSamples)];
        rel_int_struct = InterferenceFullStruct(rel_vec_idx);
        rel_powers= vertcat(rel_int_struct.Power);
        rel_freqs = vertcat(rel_int_struct.Freqs);
        min_freq = min(rel_freqs(rel_powers > 1));
        if ~isempty(min_freq)
            InterMinFreq(m) =  min_freq;
        else
            InterMinFreq(m) = randi([0,1])*50;  % if no grequency found- random between max freq and zero
        end %  if ~isempty(min_freq)
    end % for m
else % if IsAnyPreyMasked
    InterMinFreq = 0;
end % if IsAnyPreyMasked


%% Output
if ~isempty(unmasked_linked_delays)
    unmasked_linked_time = round((unmasked_linked_delays + tx_time) / filterbank_fs_ts);
else
    unmasked_linked_time = [];
end % if unmasked_linked_delays

all_zeros_vec = zeros(size(det_prey_nomaksing));
FindsMaskingStruct = struct(...
         'PulseNum', CurrPulseNum, ... v
        'DetectedPreys', det_prey_nomaksing,...v
        'DetectedPreysRxPower', nm_prey_powerdb, ...v 
        'IsAnyPreyMasked', IsAnyPreyMasked, ...v 
        'MaskedPreys', masked_prey, ... v
        'masked_prey_idx', masked_prey_idx, ...
        'TotalMaskingTimes', TotalMaskingTimes,... v
        'IsCoverTimeJamFlag' , all_zeros_vec,...v 
        'CoverTimeJamTimes' , all_zeros_vec, ...v
        'IsFwdMaskingJamFlag', all_zeros_vec,... v
        'FwdMaskingJamTimes', all_zeros_vec ,... v
        'DetecetedPrey2InterferenceRatioDB', DetecetedPrey2InterferenceRatioDB, ... %DetecetedPrey2InterferenceRatioDB, ... XXX
        'SelfCorrealationMaxdB', all_zeros_vec, ... 
        'InterCorrelationMaxdB', all_zeros_vec, ... 
        'Ref_MaskedPreys', [],  ... 
        'Ref_TotalMaskingTimes', [], ...
        'Ref_FwdMaskingJamTimes', [], ...
        'Ref_DetecetedPrey2InterferenceRatioDB', [], ... %RefDetecetedPrey2InterferenceRatioDB,...XXX
        'DetectionsDFerror', DFErr, ... %DFErr,...XXX
        'DetectionsRangeErr', RangeErr, ... % RangeErr,...XXX
        'DetectionsRelativeDirectionErr', RelativeDirectionErr, ... %RelativeDirectionErr, ...XXX
        'InterMinFreq', InterMinFreq, ...InterMinFreq, ... XXX
        'InterBW',  0,... InterBW, ...XXX
        'InterBatNum', 0, ...
        'FB_unmasked_prey', unmasked_prey, ...
        'FB_unmasked_delays', unmasked_linked_delays, ...
        'FB_unmasked_times', unmasked_linked_time, ...
        'FB_detected_Prey_times', round(nm_prey_times/ filterbank_fs_ts), ....
        'FB_detected_masking_delays', jammer_delays, ...
        'FB_detected_masking_powers', jammer_powers, ...
        'FB_estimated_masker', estimated_tx_masking_bat ...
        ); ...InterBatNum ); ... XXX