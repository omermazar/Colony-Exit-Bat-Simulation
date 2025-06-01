function [jamm_flg, pre_detected_rxLvl, max_relevant_interference_level, max_relevant_cross_interference] = ...
    Bwd_Fwd_Masking_Test(...
    wanted_detected_times, wanted_detected_pks, ...
    fwd_ntime, bwd_ntime, fwd_TH, bwd_TH, ...
    Rx_Acoustic_WantedSig, Rx_Acoustic_Interference, ...
    Inter_lvl, fb_flag, int_pks_times)

% this function calculate the SIR of the signal and reports jammong due to
% forward and backward masking
% for Correaltion Detector - the Inter_Lvl is the correlation between the
% trasmitted call and the sum of all interference signals (without the
% wanted)
% for the FilterBank the Inter_lvl is the sum of all the chhanels, shifted
% in time. and the 
% 
% Chenged Sep2022- for FilterBankDetection - 
% Inter_lvl - is the peaks of all detections in FilterBank_RxAll that are not wanted in
% times : int_pks_times

if nargin < 10
    fb_flag = 0;
    int_pks_times =[];
end % nargin    

ndetections = numel(wanted_detected_times);
max_relevant_cross_interference = zeros(1, ndetections);
max_relevant_interference_level = zeros(1, ndetections);
pre_detected_rxLvl = zeros(1, ndetections);
jamm_flg = zeros(1, ndetections);
ntot = numel(Rx_Acoustic_WantedSig);

% Inter_lvl = envelope(Inter_lvl,64,'analytic');
nlvl = 0.05* min(wanted_detected_pks);

for iWanted  = 1:ndetections
    iwanted_time = wanted_detected_times(iWanted);
    if iwanted_time < 2
        iwanted_time = 2;
    end
    iwanted_pk = wanted_detected_pks(iWanted);
    
    min_t   =  max(2,iwanted_time - fwd_ntime);
    max_t   =  min(ntot,iwanted_time + bwd_ntime);
    
    if  ~fb_flag %%% Correlation
        ix_fwd  = min_t : iwanted_time;
        fwd_lvl = max(Inter_lvl(ix_fwd));
        ix_bwd  = (iwanted_time+1):max_t;
        bwd_lvl = max(abs(Inter_lvl(ix_bwd)));

        %%% calculatating the maximum levels of the signals for SIR and Cross-Corellation
        pre_detected_rxLvl(iWanted) = max(abs(Rx_Acoustic_WantedSig([ix_fwd, ix_bwd])) );
        max_relevant_interference_level(iWanted) = max(abs(Rx_Acoustic_Interference([ix_fwd, ix_bwd])) );
        max_relevant_cross_interference(iWanted) = max(bwd_lvl, fwd_lvl);

        %%% filter-bank
        % Changed SEp 2022
        %     if fb_flag
        %         iwanted_pk = max(wanted_sum_lvl([ix_fwd, ix_bwd]));
        %         int_inrage_flag = any(int_pks_times > min_t &  int_pks_times < max_t);
        %     end % if fb_flag

        % FWD and BWD Masking
        fwd_jam_flg = (iwanted_pk ./ fwd_lvl ) < fwd_TH ; % Logical
        bwd_jam_flg = (iwanted_pk ./ bwd_lvl ) < bwd_TH ; % Logical
        jamm_flg(iWanted) = (fwd_jam_flg | bwd_jam_flg);
    else %%% FilteBank
        % Comare the detction Histogarms
        ix_fwd  = int_pks_times >= min_t & int_pks_times <= iwanted_time;
        fwd_lvl = max(Inter_lvl(ix_fwd));
        ix_bwd  = int_pks_times > iwanted_time & int_pks_times <= max_t;     (iwanted_time+1):max_t;
        bwd_lvl = max(abs(Inter_lvl(ix_bwd)));
        if ~isempty(fwd_lvl)
            fwd_jam_flg = (iwanted_pk+nlvl < fwd_lvl ); % Logical
            jamm_flg(iWanted) = fwd_jam_flg;
        end
        if ~isempty(bwd_lvl)
            bwd_jam_flg = (iwanted_pk+nlvl < bwd_lvl ) ; % Logical
            jamm_flg(iWanted) = (jamm_flg(iWanted) | bwd_jam_flg);
        end
        %%%% NEW June2024
        max_relevant_interference_level(iWanted) = max(abs(Rx_Acoustic_Interference(min_t:max_t)) );
        pre_detected_rxLvl(iWanted) = max(abs(Rx_Acoustic_WantedSig(min_t:max_t)) );
        
    end % if ~fb
    %     if fb_flag
    %         jamm_flg(iWanted) = jamm_flg(iWanted) & int_inrage_flag;
    %     end
end % for iwanted  = 1:numel(pre_detected_targets)