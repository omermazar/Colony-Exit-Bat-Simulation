function [InterReportStrctOnLine] =  ...
    AnalyzeCrossInterfernceFromOnLine(BAT, AllParams, varargin)

% function [InterReportStrct] =  AnalyzeCrossInterfernce(BatDATA.Bat, AllParams);
%
% this function returns the interferfece report for the asked bat
% This version is for offline analysis - check detection reciever types
%
% function [InterReportStrct] =  AnalyzeCrossInterfernce(BatDATA.Bat, AllParams, ReceiverType)
%   ReciverType =   'MaxDetection' -  max envelope detection (default) - Interferred if the
%                       ratio between the maximum power of the Signal and the interference is lower than SIRth
%                   'MatchedFilter' - Simple model of match filter analysis:
%                       Interferred if the total energy ration between
%                       (signal plus correlation gain) and (Interferrer) is
%                       lower than  0dB
%                   'Freq detection' Interferred if the Interfernce is in
%                   the time and freq of the echo
%

%

%%% Parameters setup %%%
% General Params
SampleTime = AllParams.SimParams.SampleTime;
NumberOfPreys = AllParams.SimParams.TotalPreysNumber;
NumOfSamples = AllParams.SimParams.SimulationTime / SampleTime + 1;
NumOfTransmittedPulses = BAT.NumOfTimesSonarTransmits;
NoiseLevel = 10.^(AllParams.BatSonarParams.NoiseLeveldB/10);
FilterBankFlag = strcmp(AllParams.BatSonarParams.ReceiverType,'FilterBank');
% NumberOfBats = AllParams.SimParams.TotalBatsNumber;

% Detection Params
DetectionTH = AllParams.BatSonarParams.PulseDetectionTH;
MinSignal2Interferer = AllParams.BatSonarParams.Signal2InterferenceRatio; % dB
FwdMaskFlag = AllParams.BatSonarParams.MaskingFwdFlag;
FwdMaskingTime = round( AllParams.BatSonarParams.MaskingFwdTime /1000 / SampleTime); % from msec to samples
BckdMaskFlag = AllParams.BatSonarParams.MaskingBckdFlag;
BckdMaskingTime = round( AllParams.BatSonarParams.MaskingBckdTime /1000 / SampleTime); % from msec to samples
MaskingFwdSIRth = AllParams.BatSonarParams.MaskingFwdSIRth; % dB
MaskingBckdSIRth = AllParams.BatSonarParams.MaskingBckdSIRth; % dB


% Input DATA
% % % % % PulsePower = BAT.PulsePower; % db
EchosStruct = BAT.EchosFromPreyStruct;
AllInterPowerVec = 10*log10(abs(BAT.InterferenceVec));
DirectInterPowerVec = 10*log10(abs(BAT.InterferenceDirectVec)); % The vector of the Interferences from other bats' direct calls
EchosInterPowerVec = 10*log10(abs(BAT.InterferenceVec-BAT.InterferenceDirectVec + NoiseLevel)); % the vector of internterfernces from received echoes fro,other bat's calls

%%% round(1e-3/SampleTime);

%%% Init Counters %%%

NumberOfPreyFinds = 0;
DetectionTimes = zeros(1,NumOfTransmittedPulses*NumberOfPreys); % maximum size
DetectionPulses = zeros(1,NumOfTransmittedPulses);
NumberOfDirectInetfrerredFinds = 0;
DirectInetfrerredTimes = DetectionTimes;
DirectInetfrerredPulseNum = DetectionTimes;

PulseHuntedPreyVec = zeros(1,NumOfTransmittedPulses);
HuntedPreyInterferredTimes = PulseHuntedPreyVec;
HuntedPreyInterferredPulseNum = PulseHuntedPreyVec;
HuntedPreyNumInterferred = PulseHuntedPreyVec;
NumberOfPulsesTryingToHunt = 0;
NumberOfHuntedPreyInterferred = 0;

NumberOfEchosInetfrerredFinds = 0;
EchosInetfrerredTimes = DetectionTimes;
EchosInetfrerredPulseNum = DetectionTimes;
% FwdMasking
NumberOfFwdMaskInetfrerredFinds = 0;
FwdMaskInetfrerredTimes = DetectionTimes;
FwdMaskInetfrerredPulseNum = DetectionTimes;
% BckdMasking
NumberOfBckdMaskInetfrerredFinds =0;
BckdMaskInetfrerredTimes = DetectionTimes;
BckdMaskInetfrerredPulseNum = DetectionTimes;

SearchInterfernceTotal = 0;
ApproachInterfernceTotal = 0;
ObsManInterfernceTotal = 0;
SearchDetectTotal = 0;
ApproachDetectTotal = 0;
ObsManDetectTotal = 0;

ManueverTypeOfPulses = repmat({''},1,NumOfTransmittedPulses); % 'Foraging' / 'Hunting' / 'ObsMan'

PulseTimesVec= zeros(1,NumOfTransmittedPulses);

% Matrix of detection of each prey
PreysDetectionMatrix = zeros(NumberOfPreys , NumOfTransmittedPulses );
PreysInterferenceMatrix = zeros(NumberOfPreys , NumOfTransmittedPulses );

%%% main function %%%

%%
%%% OnLine : contstruct the data from the 'real-time' detections %%%
% Function definision - not include nans and infs
myStd  = @(x) std(x(~isinf(x)),'omitnan');
myMean = @(x) mean(x(~isinf(x)),'omitnan');

MaskingStruct = BAT.FindsMaskingStruct;
if isfield(BAT.FindsMaskingStruct, 'IsAnyPreyMasked')
    IsAnyMaskedCell = {BAT.FindsMaskingStruct.IsAnyPreyMasked}; % cell with 1's on the pulses with making, canitns empty fields
else % if isfield
    IsAnyMaskedCell = {};
end %if isfield

%%% Genenral data
PreyFindsVec= [BAT.PreyFindsStruct.DetecectedPreyWithOutInterference]; % vertcat(BAT.PreyFindsStruct.DetecectedPreyWithOutInterference); %XXX; % Vector of all the prey items detected by the bat before checking for ineterfernce
InterReportStrctOnLine.TotalNumberOfPulses = NumOfTransmittedPulses;
InterReportStrctOnLine.TotalNumberOfDetectionsWithoutInterference = length(PreyFindsVec); % without interfernce
InterReportStrctOnLine.DetectionTimesWithoutInterference = ...
    double([BAT.PreyFindsStruct.DetectedTimesWithOutInterfernece])*SampleTime; % The times of the detections
InterReportStrctOnLine.DetectionPulsesWithoutInterference = ...
    [BAT.PreyFindsStruct(find([BAT.PreyFindsStruct.IsAnyPreyDetected])).TransmittedPulseNum ]; % The Pulse numnbers that the detection acuured in

%%% total Masking data
InterReportStrctOnLine.TotalInterferenceTimes = ...
    double(nonzeros([MaskingStruct.TotalMaskingTimes]))'*SampleTime; % The times of the maskinf to the detctions
InterReportStrctOnLine.TotalInterferencePulseNum = zeros(1 , NumOfTransmittedPulses ); % init
InterReportStrctOnLine.TotalInterference = length(InterReportStrctOnLine.TotalInterferenceTimes);
InterReportStrctOnLine.DirectInetfrerredTimes =  double(nonzeros(unique([MaskingStruct.CoverTimeJamTimes])))'*SampleTime;

%%% total cluterr
InterReportStrctOnLine.TotalPreyCluttered = length([BAT.PreyFindsStruct.ClutteredPrey]);
InterReportStrctOnLine.TotalTimesOfCluttered = [BAT.PreyFindsStruct.Clutter_nTimes]*SampleTime;

% Power and SIR
% if ~FilterBankFlag
InterReportStrctOnLine.RxPowerOfDetectedPreys = [BAT.PreyFindsStruct.RxPowerOfDetectedPreys];
InterReportStrctOnLine.SIROfDetectedPreys = [BAT.PreyFindsStruct.SIROfDetectedPreys]; %[BAT.PreyFindsStruct.SIROfDetectedPreys]; XXX
InterReportStrctOnLine.SIROfDetectedPreys_ave = myMean(InterReportStrctOnLine.SIROfDetectedPreys);
InterReportStrctOnLine.SIROfDetectedPreys_std = myStd(InterReportStrctOnLine.SIROfDetectedPreys);
InterReportStrctOnLine.ReferenceSIRByPowerTH = [BAT.FindsMaskingStruct.Ref_DetecetedPrey2InterferenceRatioDB]; % XXX

% Estimation Errors
InterReportStrctOnLine.Range_est_errs = [BAT.PreyFindsStruct.DetectionsRangeErr];
max_err = 2*std(InterReportStrctOnLine.Range_est_errs);
idx_errs = abs(InterReportStrctOnLine.Range_est_errs - myMean(InterReportStrctOnLine.Range_est_errs)) <= max_err;
InterReportStrctOnLine.Range_est_err_ave = myMean(InterReportStrctOnLine.Range_est_errs(idx_errs));
InterReportStrctOnLine.Range_est_err_std = myStd(InterReportStrctOnLine.Range_est_errs(idx_errs));

InterReportStrctOnLine.DF_est_errs = [BAT.PreyFindsStruct.DetectionsDFerror];
InterReportStrctOnLine.DF_est_err_ave = myMean(InterReportStrctOnLine.DF_est_errs(idx_errs));
InterReportStrctOnLine.DF_est_err_std = myStd(InterReportStrctOnLine.DF_est_errs(idx_errs));

% Pulses DATA
AllDetectedPulseNum = InterReportStrctOnLine.DetectionPulsesWithoutInterference;
InterReportStrctOnLine.TxPowerOfDetected = ...
    [BAT.TransmittedPulsesStruct(AllDetectedPulseNum).PulsePower];
InterReportStrctOnLine.PulseWidthOfDetected = ...
    [BAT.TransmittedPulsesStruct(AllDetectedPulseNum).PulseWidth]*SampleTime;
InterReportStrctOnLine.IPIOfDetected = ...
    [BAT.TransmittedPulsesStruct(AllDetectedPulseNum).IPItoNextPulse]*SampleTime;


% Masking to Hunted Prey
InterReportStrctOnLine.NumberOfPulsesTryingToHunt = ...
    length(nonzeros([BAT.PreyFindsStruct.PreyNumToHunt])');
InterReportStrctOnLine.NumberOfMaskedHuntingPulses = ...
    sum([BAT.PreyFindsStruct.IsHuntedPreyMasked]);

% Direct and FwdMAsking
InterReportStrctOnLine.NumberOfDirectInetfrerredFinds = length(InterReportStrctOnLine.DirectInetfrerredTimes);
InterReportStrctOnLine.FwdMaskInetfrerredTimes =  double(nonzeros(unique([MaskingStruct.FwdMaskingJamTimes])))'*SampleTime;
InterReportStrctOnLine.NumberOfFwdMaskInetfrerredFinds = length(InterReportStrctOnLine.FwdMaskInetfrerredTimes );

% refernce power th detection
InterReportStrctOnLine.Reference_PowerTHTimes = [MaskingStruct.Ref_TotalMaskingTimes];
InterReportStrctOnLine.Reference_Total = numel(InterReportStrctOnLine.Reference_PowerTHTimes);

% Ratios
if InterReportStrctOnLine.TotalNumberOfDetectionsWithoutInterference > 0

    InterReportStrctOnLine.DirectInterferenceRatio = ...
        InterReportStrctOnLine.NumberOfDirectInetfrerredFinds ./ InterReportStrctOnLine.TotalNumberOfDetectionsWithoutInterference;
    InterReportStrctOnLine.FwdMaskInterferenceRatio = ...
        InterReportStrctOnLine.NumberOfFwdMaskInetfrerredFinds ./ InterReportStrctOnLine.TotalNumberOfDetectionsWithoutInterference;
    InterReportStrctOnLine.TotalInterferenceRatio = ...
        InterReportStrctOnLine.TotalInterference ./ InterReportStrctOnLine.TotalNumberOfDetectionsWithoutInterference;
    InterReportStrctOnLine.Reference_Ratio = ...
        InterReportStrctOnLine.Reference_Total ./ InterReportStrctOnLine.TotalNumberOfDetectionsWithoutInterference;
    InterReportStrctOnLine.TotalPreyClutteredRatio =  ...
        InterReportStrctOnLine.TotalPreyCluttered ./ ...
        (InterReportStrctOnLine.TotalNumberOfDetectionsWithoutInterference+ InterReportStrctOnLine.TotalPreyCluttered);
else % if InterReportStrctOnLine.TotalNumberOfDetections > 0
    InterReportStrctOnLine.DirectInterferenceRatio = nan;
    InterReportStrctOnLine.FwdMaskInterferenceRatio = nan;
    InterReportStrctOnLine.TotalInterferenceRatio = nan;
    InterReportStrctOnLine.Reference_Ratio = nan;
    InterReportStrctOnLine.TotalPreyClutteredRatio = nan;
end % if InterReportStrctOnLine.TotalNumberOfDetections > 0

if InterReportStrctOnLine.NumberOfPulsesTryingToHunt > 0
    InterReportStrctOnLine.InterferenceToHuntedRatio = ...
        InterReportStrctOnLine.NumberOfMaskedHuntingPulses ./ numel(InterReportStrctOnLine.DetectionPulsesWithoutInterference);

else % if InterReportStrctOnLine.NumberOfPulsesTryingToHunt > 0
    InterReportStrctOnLine.InterferenceToHuntedRatio = nan;
end % if InterReportStrctOnLine.NumberOfPulsesTryingToHunt > 0

for kk=1:length(IsAnyMaskedCell)
    if IsAnyMaskedCell{kk}
        InterReportStrctOnLine.TotalInterferencePulseNum(kk) = MaskingStruct(kk).PulseNum;
    end % if IsAnyMaskedCell{kk}
end % for kTimes
InterReportStrctOnLine.TotalInterferencePulseNum = nonzeros(InterReportStrctOnLine.TotalInterferencePulseNum)';
InterReportStrctOnLine.NumberOfCatches = length(BAT.CatchPreyTimes);

%%% Maskin Power Summary
InterVecDB =  10*log10(abs(BAT.InterferenceVec));
[InterReportStrctOnLine.MaskingLevelDB_Hist.Counts, InterReportStrctOnLine.MaskingLevelDB_Hist.Centers] = ...
    hist(InterVecDB);
InterReportStrctOnLine.MaskingLevelDB_ave = myMean(InterVecDB);
InterReportStrctOnLine.MaskingLevelDB_median = median(InterVecDB);
InterReportStrctOnLine.MaskingLevelDB_Q90 = quantile(InterVecDB, 0.9);

%%% Stage Analysis
StageCell = {BAT.ManueverCmdStruct(1:NumOfTransmittedPulses).ManueverStage};

% [uStages, ~ , StagesIndex] = unique(StageCell);
% [~, ~ , StagesIndex] = unique(StageCell);
uStages = {'Search' , 'Approach', 'Buzz', 'ObstacleManuver', 'AvoidBatMan', 'CaveExit'};
% uStages = unique({BAT.ManueverCmdStruct(1:BAT.NumOfTimesSonarTransmits).ManueverStage});

SumStagesStruct.Stage = uStages';
SumStagesStruct.TotalPulses = zeros(size(uStages))';
SumStagesStruct.TotalSequences = zeros(size(uStages))';
SumStagesStruct.TotalMaskingToHunted = zeros(size(uStages))';
SumStagesStruct.MaskingToHuntedRatio = zeros(size(uStages))';
SumStagesStruct.StagePulsesCell = cell(length(uStages),1);
for k=1:length(uStages)
    StagesIndex = strcmp(StageCell, uStages(k) ); %%%
    %     SumStagesStruct.TotalPulses(k) = sum(strcmp(StageCell, uStages(k) ));
    SumStagesStruct.TotalPulses(k) = sum(StagesIndex);
    %     SumStagesStruct.TotalMaskingToHunted(k) = ...
    %         sum([BAT.ManueverCmdStruct(StagesIndex == k).HuntedPreyMaskedFlag]);
    SumStagesStruct.TotalMaskingToHunted(k) = ...
        sum([BAT.ManueverCmdStruct(StagesIndex).HuntedPreyMaskedFlag]);
    if SumStagesStruct.TotalPulses(k) >0
        SumStagesStruct.MaskingToHuntedRatio(k) = ...
            SumStagesStruct.TotalMaskingToHunted(k) ./ SumStagesStruct.TotalPulses(k);
    else % if SumStagesStruct.TotalPulses(k) >0
        SumStagesStruct.MaskingToHuntedRatio(k) = nan;
    end %if SumStagesStruct.TotalPulses(k) >0

    %     TempVec =  find((StagesIndex == k));
    TempVec =  find((StagesIndex));
    SumStagesStruct.StagePulsesCell(k) = {TempVec'};

    %%% SEQUENCES - NEW June 2022
    minSeqLength = 3;
    allStagePulses = SumStagesStruct.StagePulsesCell{k};
    if any(allStagePulses)
        diffsIx= find(diff(allStagePulses) > 1)';
        startsIX = [1, (diffsIx+1)];
        endIX = [diffsIx, numel(allStagePulses)];

        sequenceStarts =  allStagePulses(startsIX)';
        sequenceEnds= allStagePulses(endIX)';
        SumStagesStruct.TotalSequences(k) = sum(sequenceEnds- sequenceStarts +1  >= minSeqLength);
    end % if any(allStagePulses)



end %  k=1:length(uStages)

InterReportStrctOnLine.SumStagesStruct = SumStagesStruct;
InterReportStrctOnLine.SumStagesTable = struct2table(SumStagesStruct);
% rename fields

%%% Success Rate
InterReportStrctOnLine.SuccessRateBuzz     = InterReportStrctOnLine.NumberOfCatches ./ ...
    InterReportStrctOnLine.SumStagesStruct.TotalSequences(contains(uStages,'Buzz'));
InterReportStrctOnLine.SuccessRateApproach = InterReportStrctOnLine.NumberOfCatches ./ ...
    InterReportStrctOnLine.SumStagesStruct.TotalSequences(contains(uStages,'Approach'));

%% Check
% InterReportStrctOnLine.HuntingSequencesSummary = HuntingSequencesAnalysis(BAT, AllParams);

%%% stage 'value' nTime in each
PulseTimesVec = [ EchosStruct.TransmittedPulseTime];
InterReportStrctOnLine.FullFlightStages = Pulses2Flight(StageCell, PulseTimesVec, NumOfSamples);

%% Filter Bank
%%%% Nov2021
% if strcmp(AllParams.BatSonarParams.ReceiverType, 'FilterBank')
%     InterReportStrctOnLine.fb_det_SIRdb = mean([BAT.FindsMaskingStruct.DetecetedPrey2InterferenceRatioDB]);
%     InterReportStrctOnLine.fb_masking_powerdb = mean([BAT.FindsMaskingStruct.FB_detected_masking_powers]);
%     InterReportStrctOnLine.fb_mean_err_microsec = ...
%         mean([BAT.FilterBank_all_detecetd_echoes.estimated_errors])/AllParams.BatSonarParams.FilterBank_Fs/1e-6;
%     InterReportStrctOnLine.fb_std_err_microsec = ...
%         std([BAT.FilterBank_all_detecetd_echoes.estimated_errors])/AllParams.BatSonarParams.FilterBank_Fs/1e-6;
%
%     if numel(InterReportStrctOnLine.DetectionPulsesWithoutInterference) > 0
%         InterReportStrctOnLine.FA_rate = numel(nonzeros([BAT.ManueverCmdStruct.JamFalseAlarmFlag])) ./  ...
%             numel(InterReportStrctOnLine.DetectionPulsesWithoutInterference);
%     else
%         InterReportStrctOnLine.FA_rate = 0;
%     end
% end % if strcmp(AllParams.BatSonarParams.ReceiverType, 'FilterBank')
% end % function

%%

%% parameters for hunted prey
SIR2_hunted = zeros(1,NumOfTransmittedPulses);
Range2_hunted_err = SIR2_hunted;
DF2_hunted_err = SIR2_hunted;
fb_SIR2_hunted  = SIR2_hunted;
fb_err2_hunted = SIR2_hunted;


for nn = 1:max([BAT.PreyFindsStruct.TransmittedPulseNum])
    if BAT.PreyFindsStruct(nn).PreyNumToHunt > 0
        id_hunt = BAT.PreyFindsStruct(nn).DetectedPreyNum == BAT.PreyFindsStruct(nn).PreyNumToHunt;
        if sum(id_hunt)
            SIR2_hunted(nn) =  BAT.PreyFindsStruct(nn).SIROfDetectedPreys(id_hunt);
            Range2_hunted_err(nn) = BAT.PreyFindsStruct(nn).DetectionsRangeErr(id_hunt);
            DF2_hunted_err(nn) = BAT.PreyFindsStruct(nn).DetectionsDFerror(id_hunt);

        end %  if ~isempty(id_hunt)
    end % if BAT.PreyFindsStruct(nn).PreyNumToHunt > 0
end % for nn

InterReportStrctOnLine.SIR2_hunted = SIR2_hunted;
InterReportStrctOnLine.Range2_hunted_err = Range2_hunted_err;
InterReportStrctOnLine.DF2_hunted_err = DF2_hunted_err;
InterReportStrctOnLine.fb_SIR2_hunted = fb_SIR2_hunted;
InterReportStrctOnLine.fb_err2_hunted = fb_err2_hunted;

InterReportStrctOnLine.SIR2_hunted_ave = myMean(nonzeros(SIR2_hunted));
InterReportStrctOnLine.Range2_hunted_err_ave = myMean(nonzeros(Range2_hunted_err));
InterReportStrctOnLine.Range2_hunted_err_std = myStd(nonzeros(Range2_hunted_err));
InterReportStrctOnLine.DF2_hunted_err_ave = myMean(nonzeros(DF2_hunted_err));
InterReportStrctOnLine.DF2_hunted_err_std = myStd(nonzeros(DF2_hunted_err));

%% June22 -Cave Exit and  obstacle finds
% if Detect Obstacles than Aply the statsics of masking to the ObsFinds Struct
% Overwrite the relevant DATA

if AllParams.SimParams.DetectObs
    InterReportStrctOnLine.obsTotalNumberOfDetectionsWithoutInterference = ...
        numel([BAT.ObsFindsStruct.DetecectedPreyWithOutInterference]); % without interfernce
    InterReportStrctOnLine.obsDetectionTimesWithoutInterference          = ...
        double([BAT.ObsFindsStruct.DetectedTimes])*SampleTime; % The times of the detections
    InterReportStrctOnLine.obsDetectionPulsesWithoutInterference         = ...
        [BAT.ObsFindsStruct(find([BAT.ObsFindsStruct.IsAnyPreyDetected])).TransmittedPulseNum ]; % The Pulse numnbers that the detection acuured in

    %%% total Masking data
    InterReportStrctOnLine.obsTotalInterferenceTimes    = ...
        double(nonzeros([BAT.Obs_MaskingStruct.TotalMaskingTimes]))'*SampleTime; % The times of the maskinf to the detctions
    ixPulsesWithoutMasking = cellfun(@isempty, {BAT.ObsFindsStruct.MaskedPreys});
    InterReportStrctOnLine.obsTotalInterferencePulseNum = ...
        [BAT.ObsFindsStruct(~ixPulsesWithoutMasking).TransmittedPulseNum ]; % init
    InterReportStrctOnLine.obsTotalInterference         = ...
        numel([BAT.ObsFindsStruct.MaskedPreys]);
    %     InterReportStrctOnLine.DirectInetfrerredTimes =  double(nonzeros(unique([MaskingStruct.CoverTimeJamTimes])))'*SampleTime;

    % Power and SIR
    InterReportStrctOnLine.obsRxPowerOfDetected     = [BAT.ObsFindsStruct.RxPowerOfDetectedPreys];
    InterReportStrctOnLine.obsRxPowerOfDetected_ave = myMean(InterReportStrctOnLine.obsRxPowerOfDetected);
    InterReportStrctOnLine.obsRxPowerOfDetected_std = myStd(InterReportStrctOnLine.obsRxPowerOfDetected);

    % Estimation Errors
    InterReportStrctOnLine.obsRange_est_errs = [BAT.ObsFindsStruct.DetectionsRangeErr];
    max_err = 2*std(InterReportStrctOnLine.obsRange_est_errs);
    idx_errs = abs(InterReportStrctOnLine.obsRange_est_errs - myMean(InterReportStrctOnLine.obsRange_est_errs)) <= max_err;
    InterReportStrctOnLine.obsRange_est_err_ave = myMean(InterReportStrctOnLine.obsRange_est_errs(idx_errs));
    InterReportStrctOnLine.obsRange_est_err_std = std(InterReportStrctOnLine.obsRange_est_errs(idx_errs));

    InterReportStrctOnLine.obsDF_est_errs = [BAT.ObsFindsStruct.DetectionsDFerror];
    InterReportStrctOnLine.obsDF_est_err_ave = myMean(InterReportStrctOnLine.obsDF_est_errs(idx_errs));
    InterReportStrctOnLine.obsDF_est_err_std = myStd(InterReportStrctOnLine.obsDF_est_errs(idx_errs));

    % Ratios
    if InterReportStrctOnLine.obsTotalNumberOfDetectionsWithoutInterference > 0
        InterReportStrctOnLine.obsTotalInterferenceRatio = ...
            InterReportStrctOnLine.obsTotalInterference ./ InterReportStrctOnLine.obsTotalNumberOfDetectionsWithoutInterference;
    else % if InterReportStrctOnLine.TotalNumberOfDetections > 0
        InterReportStrctOnLine.obsTotalInterferenceRatio = nan;
    end % if InterReportStrctOnLine.TotalNumberOfDetections > 0
else
    InterReportStrctOnLine.obsTotalInterferenceRatio = nan;
    InterReportStrctOnLine.obsTotalInterference      = nan;
    InterReportStrctOnLine.obsRange_est_err_ave      = nan;
    InterReportStrctOnLine.obsRange_est_err_std      = nan;
    InterReportStrctOnLine.obsDF_est_err_ave         = nan;
    InterReportStrctOnLine.obsDF_est_err_std         = nan;
    InterReportStrctOnLine.obsRxPowerOfDetected_ave  = nan;
    InterReportStrctOnLine.obsRxPowerOfDetected_std  = nan;
    InterReportStrctOnLine.obsTotalNumberOfDetectionsWithoutInterference = nan;

end % if AllParams.SimParams.DetectObs

 %% Jan2023 - Add Masking Obs Decision False Alarm
nCalls = BAT.NumOfTimesSonarTransmits;
nn = zeros(1,nCalls);
for k = 1:nCalls
    nn(k) = strcmp(BAT.ManueverCmdStruct(k).ManueverType, 'ObsManByMasking');
end
InterReportStrctOnLine.obsManueverFalseAlarm_total = sum(nn);
InterReportStrctOnLine.obsManueverFalseAlarm_ratio = sum(nn)/nCalls ;

%% Conspecifs finds
if AllParams.SimParams.DetectConsps
    InterReportStrctOnLine.conspsTotalNumberOfDetectionsWithoutInterference = ...
        numel([BAT.Consps_FindsStruct.DetecectedPreyWithOutInterference ]); % without interfernce
    InterReportStrctOnLine.conspsDetectionTimesWithoutInterference          = ...
        double([BAT.Consps_FindsStruct.DetectedTimes])*SampleTime; % The times of the detections
    InterReportStrctOnLine.conspsDetectionPulsesWithoutInterference         = ...
        [BAT.Consps_FindsStruct(find([BAT.Consps_FindsStruct.IsAnyPreyDetected])).TransmittedPulseNum ]; % The Pulse numnbers that the detection acuured in

    %%% total Masking data
    InterReportStrctOnLine.conspsTotalInterferenceTimes    = ...
        double(nonzeros([BAT.Consps_MaskingStruct.TotalMaskingTimes]))'*SampleTime; % The times of the maskinf to the detctions
    ixPulsesWithoutMasking = cellfun(@isempty, {BAT.Consps_FindsStruct.MaskedPreys});
    InterReportStrctOnLine.conspsTotalInterferencePulseNum = ...
        [BAT.Consps_FindsStruct(~ixPulsesWithoutMasking).TransmittedPulseNum ]; % init
    InterReportStrctOnLine.conspsTotalInterference         = ...
        numel([BAT.Consps_FindsStruct.MaskedPreys]);
    %     InterReportStrctOnLine.DirectInetfrerredTimes =  double(nonzeros(unique([MaskingStruct.CoverTimeJamTimes])))'*SampleTime;

    % Power and SIR
    InterReportStrctOnLine.conspsRxPowerOfDetected = [BAT.Consps_FindsStruct.RxPowerOfDetectedPreys];
    InterReportStrctOnLine.conspsRxPowerOfDetected_ave = myMean(InterReportStrctOnLine.conspsRxPowerOfDetected);
    InterReportStrctOnLine.conspsRxPowerOfDetected_std = myStd(InterReportStrctOnLine.conspsRxPowerOfDetected);

    % Estimation Errors
    InterReportStrctOnLine.conspsRange_est_errs = [BAT.Consps_FindsStruct.DetectionsRangeErr];
    max_err = 2*std(InterReportStrctOnLine.conspsRange_est_errs);
    idx_errs = abs(InterReportStrctOnLine.conspsRange_est_errs - myMean(InterReportStrctOnLine.conspsRange_est_errs)) <= max_err;
    InterReportStrctOnLine.conspsRange_est_err_ave = myMean(InterReportStrctOnLine.conspsRange_est_errs(idx_errs));
    InterReportStrctOnLine.conspsRange_est_err_std = myStd(InterReportStrctOnLine.conspsRange_est_errs(idx_errs));

    InterReportStrctOnLine.conspsDF_est_errs = [BAT.Consps_FindsStruct.DetectionsDFerror];
    InterReportStrctOnLine.conspsDF_est_err_ave = myMean(InterReportStrctOnLine.conspsDF_est_errs(idx_errs));
    InterReportStrctOnLine.conspsDF_est_err_std = myStd(InterReportStrctOnLine.conspsDF_est_errs(idx_errs));

    % Ratios
    if InterReportStrctOnLine.conspsTotalNumberOfDetectionsWithoutInterference > 0
        InterReportStrctOnLine.conspsTotalInterferenceRatio = ...
            InterReportStrctOnLine.conspsTotalInterference ./ InterReportStrctOnLine.conspsTotalNumberOfDetectionsWithoutInterference;
    else % if InterReportStrctOnLine.TotalNumberOfDetections > 0
        InterReportStrctOnLine.conspsTotalInterferenceRatio = nan;
    end % if InterReportStrctOnLine.TotalNumberOfDetections > 0
else
    InterReportStrctOnLine.conspsTotalInterferenceRatio = nan;
    InterReportStrctOnLine.conspsTotalInterference      = nan;
    InterReportStrctOnLine.conspsRange_est_err_ave      = nan;
    InterReportStrctOnLine.conspsRange_est_err_std      = nan;
    InterReportStrctOnLine.conspsDF_est_err_ave         = nan;
    InterReportStrctOnLine.conspsDF_est_err_std         = nan;
    InterReportStrctOnLine.conspsRxPowerOfDetected_ave  = nan;
    InterReportStrctOnLine.conspsRxPowerOfDetected_std  = nan;
    InterReportStrctOnLine.conspsTotalNumberOfDetectionsWithoutInterference = nan;
end % if AllParams.SimParams.Detectconsps

%% %% June22 -ManuverDATA
InterReportStrctOnLine.TotalFilghtDistance = sum(sqrt(diff(BAT.xBati).^2+diff(BAT.yBati).^2), 'omitnan');
InterReportStrctOnLine.ExitSuccess         = BAT.ExitSuccess;
InterReportStrctOnLine.ExitnTimes          = BAT.ExitnTime;
InterReportStrctOnLine.CrushesObsTotal     = BAT.CrushesObsNum;
InterReportStrctOnLine.CrushesObsnTimes    = BAT.CrushesObsnTimes;
InterReportStrctOnLine.CrushesConspsTotal  = BAT.CrushesConspsNum;
InterReportStrctOnLine.CrushesConspsnTimes = BAT.CrushesConspsnTimes;

% General Flight DAta
InterReportStrctOnLine.xBat                   = BAT.xBatPos;
InterReportStrctOnLine.yBat                   = BAT.yBatPos;
InterReportStrctOnLine.velocity               = BAT.BatVelocity / AllParams.SimParams.SampleTime * AllParams.SimParams.xyResolution;
InterReportStrctOnLine.dist2Consps_Mat        = reshape([BAT.ConspLocStruct.Distances]*AllParams.SimParams.xyResolution, AllParams.SimParams.TotalBatsNumber-1,[]);
InterReportStrctOnLine.dist2Consps_MeanByPulse = mean(InterReportStrctOnLine.dist2Consps_Mat, 'omitnan');
InterReportStrctOnLine.dist2Consps_STDByPulse  = std(InterReportStrctOnLine.dist2Consps_Mat, 'omitnan');
InterReportStrctOnLine.dist2Consps_MinByPulse  = min(InterReportStrctOnLine.dist2Consps_Mat);
InterReportStrctOnLine.dist2Consps_MedByPulse  = median(InterReportStrctOnLine.dist2Consps_Mat, 'omitnan');
InterReportStrctOnLine.dist2Consps_MeanMin     = myMean(InterReportStrctOnLine.dist2Consps_MinByPulse);
InterReportStrctOnLine.dist2Consps_STDMin      = myStd(InterReportStrctOnLine.dist2Consps_MinByPulse);
InterReportStrctOnLine.dist2Consps_MeanMean    = myMean(InterReportStrctOnLine.dist2Consps_MeanByPulse);
InterReportStrctOnLine.dist2Consps_STDMean     = myStd(InterReportStrctOnLine.dist2Consps_MeanByPulse);
InterReportStrctOnLine.dist2Consps_MeanMed     = myMean(InterReportStrctOnLine.dist2Consps_MedByPulse);
InterReportStrctOnLine.dist2Consps_STDMed      = myStd(InterReportStrctOnLine.dist2Consps_MedByPulse);


% Releveant Echoloaction DAta
InterReportStrctOnLine.EchosFromPreyStruct  = BAT.EchosFromPreyStruct;
if AllParams.SimParams.DetectObs
    InterReportStrctOnLine.EchosFromObsStruct   = BAT.EchosFromObsStruct;
end
if AllParams.SimParams.DetectConsps
    InterReportStrctOnLine.EchosFromConspStruct = BAT.EchosFromConspStruct;
end
InterReportStrctOnLine.TransmittedPulses    = BAT.TransmittedPulsesStruct;
InterReportStrctOnLine.PreyFindsStruct      = BAT.PreyFindsStruct;
InterReportStrctOnLine.Consps_FindsStruct   = BAT.Consps_FindsStruct;
InterReportStrctOnLine.ObsFindsStruct       = BAT.ObsFindsStruct;

%% Probabilty of detection in Distancece to Abppraoch between +-60 detgrees infront of the bat

% detDist  = AllParams.BatFlightParams.DistanceFromPreyToApproach;
detDistTest  = [1, 2, 3]; % meters
detAngle = pi/3;

for ii = 1:numel(detDistTest)
    detDist = detDistTest(ii);

    %%% Prey %%%
    if AllParams.SimParams.DetectPrey
        DetAllDist  = [BAT.PreyFindsStruct.Dist2DetectedPrey] * AllParams.SimParams.xyResolution;
        DetAllAngle = [BAT.PreyFindsStruct.Angle2DetectedPrey];
        % the distances
        allDist  = [BAT.Vector2Preys.Dist2Prey];
        allAngle = [BAT.Vector2Preys.Angle2Prey];
        ixCalls  = diff(allDist) ~= 0;
        allDist  = allDist(ixCalls);
        allAngle = allAngle(ixCalls);

        ixDetDet = DetAllDist < detDist & abs(DetAllAngle) < detAngle;
        ixDetAll = allDist < detDist & abs(allAngle) < detAngle;

        InterReportStrctOnLine.(['preyTotalInBeam_', num2str(detDist),'m']) = sum(ixDetAll);
        if sum(ixDetAll) > 0
            InterReportStrctOnLine.(['preyDetectProb_', num2str(detDist),'m']) = sum(ixDetDet) ./ sum(ixDetAll);
        else
            InterReportStrctOnLine.(['preyDetectProb_', num2str(detDist),'m']) = nan;
        end % if sum(ixDetAll) > 0

    else
        InterReportStrctOnLine.(['preyDetectProb_', num2str(detDist),'m'])  = nan;
        InterReportStrctOnLine.(['preyTotalInBeam_', num2str(detDist),'m']) = nan;

    end % if AllParams.SimParams.DetectPrey


    %%% Clutter %%%
    if AllParams.SimParams.DetectObs
        DetAllDist  = [BAT.ObsFindsStruct.Dist2DetectedPrey] * AllParams.SimParams.xyResolution;
        DetAllAngle = [BAT.ObsFindsStruct.Angle2DetectedPrey];
        % the distances
        allDist  = [BAT.ObsInBeamStruct.Distances] * AllParams.SimParams.xyResolution;
        allAngle = [BAT.ObsInBeamStruct.TargetAngle];

        ixDetDet = DetAllDist < detDist & abs(DetAllAngle) < detAngle;
        ixDetAll = allDist < detDist & abs(allAngle) < detAngle;
        % detectProb = numel(detAllDist) ./ sum(ixDetAll);
        InterReportStrctOnLine.(['obsTotalInBeam_', num2str(detDist),'m']) = sum(ixDetAll);
        if sum(ixDetAll) > 0
            InterReportStrctOnLine.(['obsDetectProb_', num2str(detDist),'m']) = sum(ixDetDet) ./ sum(ixDetAll);
        else
            InterReportStrctOnLine.(['obsDetectProb_', num2str(detDist),'m']) = nan;
        end % if sum(ixDetAll) > 0

    else % if AllParams.SimParams.DetectObs
        InterReportStrctOnLine.(['obsDetectProb_', num2str(detDist),'m'])  = nan;
        InterReportStrctOnLine.(['obsTotalInBeam_', num2str(detDist),'m']) = nan;
    end % if AllParams.SimParams.DetectObs

    %%% Consps %%%
    if AllParams.SimParams.DetectConsps
        DetAllDist  = [BAT.Consps_FindsStruct.Dist2DetectedPrey] * AllParams.SimParams.xyResolution;
        DetAllAngle = [BAT.Consps_FindsStruct.Angle2DetectedPrey];
        % the distances
        allDist  = [BAT.OtherBatsPolarCoordinates.Distances] * AllParams.SimParams.xyResolution;
        allAngle = [BAT.OtherBatsPolarCoordinates.TargetAngle];

        ixDetDet = DetAllDist < detDist & abs(DetAllAngle) < detAngle;
        ixDetAll = allDist < detDist & abs(allAngle) < detAngle;
        % detectProb = numel(detAllDist) ./ sum(ixDetAll);
        InterReportStrctOnLine.(['conspsTotalInBeam_', num2str(detDist),'m']) = sum(ixDetAll);
        if sum(ixDetAll) > 0
            InterReportStrctOnLine.(['conspsDetectProb_', num2str(detDist),'m']) = sum(ixDetDet) ./ sum(ixDetAll);
        else
            InterReportStrctOnLine.(['conspsDetectProb_', num2str(detDist),'m']) = nan;
        end % if sum(ixDetAll) > 0

    else % if AllParams.SimParams.DetectConsps
        InterReportStrctOnLine.(['conspsDetectProb_', num2str(detDist),'m'])  = nan;
        InterReportStrctOnLine.(['conspsTotalInBeam_', num2str(detDist),'m'])= nan;

    end % if AllParams.SimParams.DetectConsps

end % for ii = 1:numel(detDistAll)

