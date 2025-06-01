function [filter_bank_strct] = FilterBank_BuildStruct(BAT, CurrentPulse, FilterBank, fs_ts)

% calculates the time and freqs of the tx signal and rx echoes in high sampling rate 
% 
%   struct with fields: (example)
% 
%            PulseNum: 1
%            PulsePower: 1.0000e+11
%            TxTimes_us: [1×1750 double]
%           TxStartTime: 125
%         PulseDuration: 1750
%          TxFrequncies: [1×1750 double]
%           TxActiveFcs: [1×1750 double]
%       TxActiveFcs_idx: [1750×1 double]
%     Tx_Fcs_StartTimes: [125 250 375 500 625 750 875 1000 1125 1250 1375 1500 1625 1750]
%       RxEchosFromPrey: [1×2 struct]
%   fs_ts =  AllParams.BatSonarParams.FilterBank_Fs * AllParams.SimParams.SampleTime;

DetailedEchosStrct = BAT.EchosFromPreyStruct(CurrentPulse).EchoDetailed; 
BAT.CurrentPulseDuration =  BAT.TransmittedPulsesStruct(CurrentPulse).PulseDuration;

filter_bank_strct.PulseNum = CurrentPulse;
filter_bank_strct.PulsePower = 10.^(BAT.TransmittedPulsesStruct(CurrentPulse).PulsePower/10);
filter_bank_strct.TxTimes_us = ...
                       BAT.TransmittedPulsesStruct(CurrentPulse).StartPulseTime * fs_ts : ...
                       (BAT.TransmittedPulsesStruct(CurrentPulse).StartPulseTime+  BAT.TransmittedPulsesStruct(CurrentPulse).PulseDuration)*fs_ts-1;

filter_bank_strct.TxStartTime = filter_bank_strct.TxTimes_us(1);
filter_bank_strct.PulseDuration = BAT.CurrentPulseDuration* fs_ts;
filter_bank_strct.PulseIPI = BAT.TransmittedPulsesStruct(CurrentPulse).IPItoNextPulse * fs_ts;

Tx_freq_times = linspace( filter_bank_strct.TxStartTime, filter_bank_strct.TxTimes_us(end),...
    max(2, BAT.CurrentPulseDuration) );  
Tx_freqs =  BAT.TransmittedPulsesStruct(CurrentPulse).PulseFreqsCommands*1000; %Hz % XXXXX
if numel(Tx_freqs) < 2
    Tx_freqs =  BAT.TransmittedPulsesStruct(CurrentPulse).ChirpMinMaxFreqs([2,1])*1000; %Hz
end % if numel
filter_bank_strct.TxFrequncies = interp1(Tx_freq_times, Tx_freqs, filter_bank_strct.TxTimes_us);

[fcs_Index,~] = dsearchn(FilterBank.filter_fc', filter_bank_strct.TxFrequncies');
filter_bank_strct.TxActiveFcs = FilterBank.filter_fc(fcs_Index);
filter_bank_strct.TxActiveFcs_idx = fcs_Index;
filter_bank_strct.Tx_Fcs_StartTimes =  [BAT.TransmittedPulsesStruct(CurrentPulse).StartPulseTime : ...
                       (BAT.TransmittedPulsesStruct(CurrentPulse).StartPulseTime+  BAT.CurrentPulseDuration-1)]*fs_ts;
                                      
filter_bank_strct.RxEchosFromPrey(max(numel(DetailedEchosStrct),1)) = struct( ...
    'target_idx',[], ...
    'Rx_powers_lin',[], ...
    'Rx_times',[], ...
    'Rx_freqs',[]);

if numel(DetailedEchosStrct) >= 1
    for k_echo = 1:numel(DetailedEchosStrct)
         filter_bank_strct.RxEchosFromPrey(k_echo).target_idx = DetailedEchosStrct(k_echo).TargetIndex;
% %          filter_bank_strct.RxEchosFromPrey(k_echo).Rx_times = filter_bank_strct.TxStartTime + ...
% %              round(BAT.EchosFromPreyStruct(CurrentPulse).EchosTimes(k_echo) + ([0:(BAT.CurrentPulseDuration-1)])* fs_ts);
            
         echo_start_time = filter_bank_strct.TxStartTime + BAT.EchosFromPreyStruct(CurrentPulse).EchosTimes(k_echo)*fs_ts;
         RxTimes_ds = linspace( echo_start_time, echo_start_time + filter_bank_strct.PulseDuration, ...
              max(2, BAT.CurrentPulseDuration));
         
         filter_bank_strct.RxEchosFromPrey(k_echo).Rx_times = echo_start_time: (echo_start_time + BAT.CurrentPulseDuration*fs_ts-1); 
         Rx_powers_lin_ds = DetailedEchosStrct(k_echo).EchoAttenuationFull * filter_bank_strct.PulsePower;
         if numel(Rx_powers_lin_ds) ==1
             Rx_powers_lin_ds = [Rx_powers_lin_ds, Rx_powers_lin_ds];
         end %if numel
         filter_bank_strct.RxEchosFromPrey(k_echo).Rx_powers_lin = interp1(...
             RxTimes_ds, Rx_powers_lin_ds,...
             filter_bank_strct.RxEchosFromPrey(k_echo).Rx_times);
% %         [us_zoh_mat] = up_sampling_zoh(numel(RxTimes_ds), fs_ts);
% %         filter_bank_strct.RxEchosFromPrey(k_echo).Rx_powers_lin = Rx_powers_lin_ds * us_zoh_mat; 
        
         filter_bank_strct.RxEchosFromPrey(k_echo).Rx_freqs = filter_bank_strct.TxFrequncies;    
% %         filter_bank_strct.RxEchosFromPrey(k_echo).target_idx = DetailedEchosStrct(k_echo).TargetIndex;
% %         filter_bank_strct.RxEchosFromPrey(k_echo).Rx_powers_lin = DetailedEchosStrct(k_echo).EchoAttenuationFull * filter_bank_strct.PulsePower;
% %         filter_bank_strct.RxEchosFromPrey(k_echo).Rx_times = filter_bank_strct.TxStartTime + ...
% %             round((BAT.EchosFromPreyStruct(CurrentPulse).EchosTimes(k_echo) + [0:(BAT.CurrentPulseDuration-1)])* fs_ts);
% %         filter_bank_strct.RxEchosFromPrey(k_echo).Rx_freqs = round(DetailedEchosStrct(k_echo).EchosFreqsFull); % change to closest ...
    end % for k_echo = 1:numel(DetailedEchosStrct)
end %if numel(DetailedEchosStrct) > =1