%%
figure
sgtitle(BatNum)
tiledlayout('flow');
nexttile; hold on
% the Signals
plot(20*log10(abs(BAT(BatNum).AcousticSig_Clutter)+1), 'g', 'DisplayName', 'Obs')
plot(20*log10(abs(BAT(BatNum).AcousticSig_Calls)+1), 'k', 'DisplayName', 'Calls')
plot(20*log10(abs(BAT(BatNum).AcousticSig_ConspsObsEchoes)+1), 'm', 'DisplayName', 'ConspObs')
plot(20*log10(abs(BAT(BatNum).AcousticSig_ConspsCalls)+1), 'r', 'DisplayName', 'ConspCall')
%%% the first second
xlim([0, FsAcoustic])

%%
plot(relevant_idx, smooth(20*log10(abs(Rx_Acoustic_WantedSig)+1)), 'b', 'DisplayName', 'Wanted')


%%% the table
plot(echoTableObs.time*FsAcoustic, 10*log10(abs(echoTableObs.receivedLevel)+1), 'ok')

%%% THe Detectrtion
%%%%  without Confusiot
% nexttile; hold on
ObFindsNo = BAT(BatNum).ObsFindsStruct_NoConfusion(BAT(BatNum).CurrPulseNum);
ixDet = ismember(ObFindsNo.DetecectedPreyWithOutInterference, ObFindsNo.DetectedPreyNum);
plot(ObFindsNo.DetectedTimesWithOutInterfernece(ixDet)*SampleTime*FsAcoustic,  ObFindsNo.RxPowerOfDetectedPreys(ixDet), '*k', "MarkerSize", 8, 'LineWidth', 1.5, 'DisplayName', 'Detections No Confusion')
ObFinds = BAT(BatNum).ObsFindsStruct(BAT(BatNum).CurrPulseNum);
ixDet = ismember(ObFinds.DetecectedPreyWithOutInterference, ObFinds.DetectedPreyNum);
plot(ObFinds.DetectedTimesWithOutInterfernece(ixDet)*SampleTime*FsAcoustic,  ObFinds.RxPowerOfDetectedPreys(ixDet), 'or', "MarkerSize", 8, 'LineWidth', 1.5, 'DisplayName', 'Detections Confusion')