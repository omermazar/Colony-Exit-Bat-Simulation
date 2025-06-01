function [ cors, times, pks, pks_times  ] = my_corr_receiver(Rx_Acoustic, Tx_Acoustic_Call, min_detect_res, corrTH, detectMode)

%  Calculates the correlations, and fins the peaks of the envelope, with several TH
% [ pks, pks_times ] = my_corr_receiver(Rx_Acoustic, Tx_Acoustic_Call, min_detect_res, corrTH, time_tol)
% return the peaks Values and Times, of the Correlation betrweeen the
% Recieved signal (Rx_Acoustic) and the transmitted call(Tx_Acoustic_Call).
% min_detect_res is the sample difference between peaks, corrTH is the
% minimum value of the peak
%%% NOV2021

%% the correlation
[cors, times] = xcorr(Rx_Acoustic, Tx_Acoustic_Call);
% Normalization
cors = 2* 1/numel(Tx_Acoustic_Call) * cors;
% Take postive time s of correation
idx0 = numel(Rx_Acoustic);
cors = cors(idx0:end);
times = times(idx0:end);

%% the Receiver
switch detectMode
    case 'envelope'
        [pks, pks_ix, ~, ~] = findpeaks(envelope(cors, 64,'analytic'), 'MinPeakDistance', min_detect_res, 'MinPeakHeight', corrTH);
    case 'none'
        [pks, pks_ix, ~, ~] = findpeaks(cors, 'MinPeakDistance', min_detect_res, 'MinPeakHeight', corrTH);
end % switch

pks_times = times(pks_ix);
