        % CASE FILTERBANK

        time_tol = 0.2*1e-3* FsAcoustic; % resolution of peaks detection for comapring
  
%         fc  = FilterBank.filter_fc; % the center freqs of the filterbank 
        %%% LPF
        f_lowpass= AllParams.BatSonarParams.FilterBank_LPF_Freq*1000;
%         f_lowpass= 40*1000;
        [LPF_filt.b_filt, LPF_filt.a_filt] = butter(6, f_lowpass/(FsAcoustic/2));
        
        %%% Calculate the refence delays of the Call 
        FilterBank_RefCall = FilterBank_Acous_Output(FilterBank, LPF_filt, Tx_Acoustic_Call, FsAcoustic, true);
        
        %%% Wanted echoes(i.e. Prey echoes) %%%%
        FilterBank_RxWanted = FilterBank_Acous_Output(FilterBank, LPF_filt, Rx_Acoustic_WantedSig, FsAcoustic, false, FilterBank_RefCall);
        wanted_pks = FilterBank_RxWanted.det_pks	;
        wanted_pks_times = FilterBank_RxWanted.det_ipks	;
        corsW = FilterBank_RxWanted.Sum_out;
        
% %         [ corsW, timesW, wanted_pks, wanted_pks_times ] = my_corr_receiver(Rx_Acoustic_WantedSig, Tx_Acoustic_Call, min_detect_res, corrTH);
        
        %%% Comparing the wanted correlation to the pre-detections
        [wanted_detected_times, wanted_detected_pks, detected_ix, time_diff_wanted] = ...
            detections_compare(pre_detect_times, wanted_pks_times, wanted_pks, time_tol);
        pre_detected_targets = DetectedPreysVec(detected_ix);
        pre_detected_rxLvldB = CurrentPick(detected_ix);
        pre_undetected_targets = DetectedPreysVec(~detected_ix);
        
         
        %%% All Rx Signals %%%%
        FilterBank_RxAll = FilterBank_Acous_Output(FilterBank, LPF_filt, Rx_Acoustic_Sig, FsAcoustic, false, FilterBank_RefCall);
        all_pks = FilterBank_RxAll.det_pks;
        all_pks_times = FilterBank_RxAll.det_ipks;
        corsA = FilterBank_RxAll.Sum_out;
        
        %%% Interference Only %%%
        FilterBank_RxInt = FilterBank_Acous_Output(FilterBank, LPF_filt, Rx_Acoustic_Interference, FsAcoustic, false, FilterBank_RefCall);
        int_pks = FilterBank_RxInt.det_pks;
        int_pks_times = FilterBank_RxInt.det_ipks;
        int_cors = FilterBank_RxInt.Sum_out;
 
        Interference_corr_lvl = envelope(int_cors,64,'analytic');
        
        %%% Direct Jamming: Comparing All Detections and Wanted
        
        [all_detected_times, all_detected_pks, all_ix, time_diff_all] = ...
            detections_compare(wanted_detected_times, all_pks_times, all_pks, time_tol);         
        jamm_flg = ~all_ix;
        
         %% Test for BWD and FWD Masking An calcute SIR
        
        [fbwd_jamm_flg, max_relevant_interference_level, max_relevant_cross_interference] = ...
            Bwd_Fwd_Masking_Test(...
            wanted_detected_times, wanted_detected_pks, ...
            FwdMasking_AcousTime, BwdMasking_AcousTime, MaskingFwd_AcousTH, MaskingBwd_AcousTH, ...
            Rx_Acoustic_WantedSig, Rx_Acoustic_Interference, ...
            int_cors);
        
        jamm_flg = jamm_flg | fbwd_jamm_flg;
        
        %%% Output 
        jammed_targets = pre_detected_targets(jamm_flg) ;
        jammed_times = wanted_detected_times(jamm_flg);
        
        IsPreyMaskedFlag = jamm_flg;
        EchoSelfCorrealationMaxdB = 20*log10(wanted_detected_pks);
        
        Corr2MaskRatio = 20*log10(wanted_detected_pks./max_relevant_cross_interference);
        EchoMaxInterCorrelationdB =  20*log10(max_relevant_cross_interference);
        
        DetecetedPrey2InterferenceRatioDB = pre_detected_rxLvldB - 10*log10(abs(max_relevant_interference_level));
        
        MaskingTimes_Acoustics = TxnTime*FsAcoustic*SampleTime + jammed_times;
        TotalMaskingTimes = round(MaskingTimes_Acoustics/(FsAcoustic*SampleTime));