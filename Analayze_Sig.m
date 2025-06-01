function [] = Analayze_Sig(BatDATA, PulseNum, BatNum, ReceiverType, DetectedObj)

% Changed Sep2022
AllParams = BatDATA.AllParams;
if nargin <5
    DetectedObj = 'Prey'; % 'Prey'; 'Obs'; 'Conspecific'
end

if nargin == 4
    AllParams.BatSonarParams.ReceiverType = ReceiverType;
end

if nargin < 3
    BatNum = 1;
end

FsAcoustic = AllParams.SimParams.FsAcoustic;
SampleTime = AllParams.SimParams.SampleTime;
AcousticNumel = FsAcoustic*AllParams.SimParams.SimulationTime + 1;
BAT = BatDATA.BAT(BatNum);
PulsePower = BAT.TransmittedPulsesStruct(PulseNum).PulsePower;

%%%% NOTICE- For now Analysisi only one of the mode - TO BE FIXED %%%%%
% 01Aug22 , 
% 08Sep22

relevant_idx = round( BAT.TransmittedPulsesStruct(PulseNum).StartPulseTime* FsAcoustic*SampleTime + ...
                                        (1: BAT.TransmittedPulsesStruct(PulseNum).IPItoNextPulse* FsAcoustic*SampleTime ) ); 
try
    switch DetectedObj
        case 'Prey'
            CurrEchosStruct   = BAT.EchosFromPreyStruct(PulseNum);
            DetectedEchoesVec = BAT.PreyFindsStruct(PulseNum).DetecectedPreyWithOutInterference;
            Rx_Acoustic_WantedSig = BAT.AcousticSig_PreyEchoes(relevant_idx) +  BAT.RandNoise(relevant_idx);
        case 'Obs'
            CurrEchosStruct   = BAT.EchosFromObsStruct(PulseNum);
            DetectedEchoesVec = BAT.ObsInBeamStruct(PulseNum).DetectetedTargetsVec;
            Rx_Acoustic_WantedSig = BAT.AcousticSig_Clutter(relevant_idx) +  BAT.RandNoise(relevant_idx);
        case 'Conspecific'
            CurrEchosStruct   = BAT.EchosFromConspStruct(PulseNum);
            DetectedEchoesVec = BAT.Consps_FindsStruct(PulseNum).DetecectedPreyWithOutInterference;
            Rx_Acoustic_WantedSig = BAT.AcousticSig_CospsEchoes(relevant_idx) +  BAT.RandNoise(relevant_idx);
    end % switch
catch
    warning('Detected Object does not exist')
    return
end % try

%%%%

if relevant_idx(end) > AcousticNumel 
    relevant_idx = relevant_idx(1):AcousticNumel;
end % if
Tx_Acoustic_Call = BAT.TransmittedPulsesStruct(PulseNum).TxAcousticSig; % the call of the bat
Rx_Acoustic_Sig = BAT.AcousticSig_All(relevant_idx); % the whole Signal including interferences and noise level

% Rx_Acoustic_WantedSig =  BAT.AcousticSig_Wanted(relevant_idx);  %  BAT.AcousticSig_PreyEchoes(relevant_idx);% the wanted Signal
% BatNumRx=1;
FilterBank = BatDATA.FilterBank;

allBATs = BatDATA.BAT;

[strct] = Interfernce2DetectedPreys_Acoustics(PulsePower,...
                        DetectedEchoesVec , CurrEchosStruct, Tx_Acoustic_Call, ...
                        Rx_Acoustic_WantedSig,  Rx_Acoustic_Sig, PulseNum, allBATs, BatNum, AllParams, FilterBank , true, DetectedObj);
%% Debug version
% strctXX = Interfernce2DetectedPreys_Acoustics_ver08Feb(PulsePower,...
%                         DetectedEchoesVec , CurrEchosStruct, Tx_Acoustic_Call, ...
%                         Rx_Acoustic_WantedSig,  Rx_Acoustic_Sig, PulseNum, allBATs, BatNum, AllParams, FilterBank , false, DetectedObj);
%%
colorbar('off')
end