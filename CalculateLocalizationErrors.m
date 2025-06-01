function [DFErr, RangeErr, RelativeDirectionErr] = ...
    CalculateLocalizationErrors(...
        BAT, CurrEchosFromPreyStruct, DetectedPreysVec, ...
        SIRDB, RxPowerDB, PulseNum, AllParams)

 %%% Calculate lcalization error (DF, Distance, RElative pos)
% CalculateLocalizationErros function nedees; the Rxpower, SIR, Angle,
% Distance, TransmittedPulseTime
% Changed 19Jan21 - Swarm detection 

%%% General parmaters
SNRdb= min([SIRDB; RxPowerDB]); % signal to noise is the min between sNR and SIR
SNR = 10.^(SNRdb/10); % linear
PulseTimeWidth = BAT.TransmittedPulsesStruct(PulseNum).PulseWidth .* AllParams.SimParams.SampleTime;

xyResolution = AllParams.SimParams.xyResolution;
SampleTime = AllParams.SimParams.SampleTime;
SoundV0 =  AllParams.SimParams.SoundV0/xyResolution*SampleTime;% Sound Velocit in xyUnits

% try % add to AllParams
    IntegrationTime = AllParams.BatSonarParams.CorrelationTimeResolution * AllParams.SimParams.SampleTime;
% catch
%    IntegrationTime = 0.1*1e-3; % 100 micro
% end% try

NumberOfMeasures = PulseTimeWidth./ IntegrationTime;

%%% New For swarm- detection - Omer Jan 21
% EchosAngles = zeros(1, max(DetectedPreysVec));
% EchosAngles(DetectedPreysVec) = CurrEchosFromPreyStruct.EchosAngles;
% EchosmDistances = zeros(1, max(DetectedPreysVec));
% EchosmDistances(DetectedPreysVec) = CurrEchosFromPreyStruct.EchosmDistances;
% 
% DetectedTargetsAngles = EchosAngles(DetectedPreysVec);
% DetectedTargetsDistances = EchosmDistances(DetectedPreysVec);
%%%

%%% Old
[~,DetectedPreysVec] = ismember(DetectedPreysVec, CurrEchosFromPreyStruct.TargetIndex); % New Nov22
DetectedTargetsAngles = CurrEchosFromPreyStruct.EchosAngles(DetectedPreysVec);
DetectedTargetsDistances = CurrEchosFromPreyStruct.EchosmDistances(DetectedPreysVec);
%%%

xyResolution = AllParams.SimParams.xyResolution;
SampleTime = AllParams.SimParams.SampleTime;
SoundV0 =  AllParams.SimParams.SoundV0/xyResolution*SampleTime;% Sound Velocit in xyUnits

%% DF error %%
% Anoise = d/lambda, d= head widfth~2cm, lamda@40kHz = 8mm
Anoise = 0.5*pi/180;% *10; % 1deg error at SNR =10db
% DFNoiseErrVec =  min (Anoise / sqrt(NumberOfMeasures) ./ SNR , pi/2 ) ; % the max error is pi/2
DFNoiseErrVec =  min (Anoise ./ SNR , pi/2 ) ; % the max error is pi/2

% the resolution of the measurement as function of teta

% blur = 0.5deg @angle= 0deg, ~10deg~angle= 90deg
a = 1*pi/180; 
b = 0.2*pi/180; % 1.5 old
DFBlurrVec = b + a.*abs(sin(DetectedTargetsAngles));

DFErrSTD = sqrt(DFNoiseErrVec.^2 + DFBlurrVec.^2);
DFErr = DFErrSTD .* randn(size(DetectedPreysVec));

 %% New Max Error Allowed 2024
 maxDFErr = 3/180*pi;
 ixTrim = find(abs(DFErr) > maxDFErr);
 if any(ixTrim)
     DFErr(ixTrim) = sign(DFErr(ixTrim)).*maxDFErr;
 end

%% Range Error %%
% error = ~ m/(bw*duration) *  1/snr, err@10dBsnr = ~2cm
% estimation Range error = M*1/snr;
K1 = 10*0.01 / xyResolution ; %1cm@SNR = 10DB 
TimeRes = 50*1e-6 ./ SampleTime; % 50 micro-sec. ref Hearing By Bats chapter 3.1
discRes = 0.5* TimeRes* SoundV0; % time measure resolution [~0.8cm] 0.5 because it is  2way ; 

RangeErrSTD = min ( sqrt((K1 ./ SNR).^2 + discRes.^2) , 100); % the max err is 1 m
RangeErr = RangeErrSTD .* randn(size(DetectedPreysVec));

% RelativeDirectionErr
% the RelativeDirectionErr is complex estimation, depend on SNR and angle
% RelativeDirectionErr = G* DFErr
G = 4;
RelativeDirectionErrSTD = min(G .*DFErrSTD , pi);
RelativeDirectionErr = RelativeDirectionErrSTD .* rand(size(DetectedPreysVec));


    