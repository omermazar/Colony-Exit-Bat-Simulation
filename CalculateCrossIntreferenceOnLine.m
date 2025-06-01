function [DirectInterfernceTimes, DirectInterferncePower, AllInetfernceTimes,...
    AllInetferncePower, ALLInterStruct, TransmitFreqs, BiStatDetectionSrtct, PhantomEchoesStrct, FilterBank_jam_struct, AcousticSig_ConspsCalls] = ...
    CalculateCrossIntreferenceOnLine...
    (BATSstrct, BATTxNum, BatRxNum, PREY, CoordinationFlag, AllParams, CurrPulseNum, FilterBank, jam_mat, varargin)
%%
%%
% this function calcultes and returns the interfence that BATTrantmitter causes to
% BATReceiver
% using the equation : Pr= Pt*Dt*Dr*(lambda/(4*pi*R))^2 * exp(alpha*R)
% where"    Pr= Pweor recieved
%           Pt= tramsmitted Power
%           Dt- Directivity of transmitter (antenna). , function pf angle
%           Dr - Directivity of reciever (antenna), function pf angle
%           lambda = wavelength = SoundVelocity/ Frequency
%           alpha- Atmosphere absorbtion parameter
%           R - distance
%           CurrPulseNum - the pulse we calculate
% OutPuts:  
%           Bat2RecievedPulseTimes - the nTimes vector of the recieved interference
%           Bat2Attenuation - the vector of the received powers 
%           EchosIntRecievedPower, EchosIntRecievedTimes - the power and
%               times of the echos of targets from the interferer bat
%           Freqs - vector of the frequncy recieved
%           BiStatDetectionSrtct - struct of potential detection form other
%               bats echoes
%

% 15/04/19 updated by omer : Adding bi-stat detection output

%%% 7/11/21 
%%% The Building of the Acoustic Signal of the Interference is
%%% executerd for each kRxBat in the function CalculateCrossIntreferenceOnLine
%%%%% NOV2021

%%
%%%% Inputs

BATTrantmitter = BATSstrct(BATTxNum);
BATReciever =  BATSstrct(BatRxNum);
%%% CONSTS and Parameters from Inputs
mSoundVelocity = AllParams.SimParams.SoundV0; %m/sec
SampleTime = AllParams.SimParams.SampleTime;
xyResolution = AllParams.SimParams.xyResolution;
MaxTime = AllParams.SimParams.SimulationTime ./ SampleTime ;
% TranmittedPulsePower = 10^(BATTrantmitter.PulsePower/10);
AttenuationConst = 1; % Just parameter
% Nov2021
FilterBank_jam_struct = struct();

%%                   
%%% the tranmitted time, power and freqs
TransmitStartTime = BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).StartPulseTime;
Bat1TransmitPulseTimes = TransmitStartTime:(TransmitStartTime+ BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).PulseDuration-1); 
TransmitFreqs = BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).PulseFreqsCommands;
TranmittedPulsePower = BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).PulsePower;
PulseDuration = length(TransmitFreqs);
%%

AlphaAtmAtt = 0.038*TransmitFreqs - 0.3; %  from table, @20Cel.50% humidity Atmosphere attenuation[m^(-1)] ...
                    %  1dB loose/meter % 0.3 for 40kHa (at 30C,70% humidty],alpha =  0.1 for 28kH
%%
%%% the position of the bats when Bat1 transmits
xBatTx = deal(BATTrantmitter.xBati(TransmitStartTime));
yBatTx = deal(BATTrantmitter.yBati(TransmitStartTime));
BatTxDirection = deal(BATTrantmitter.Teta(TransmitStartTime));
xBatRx = deal(BATReciever.xBati(TransmitStartTime));
yBatRx = deal(BATReciever.yBati(TransmitStartTime));
BatRxDirection = deal(BATReciever.Teta(TransmitStartTime));

%%% the distances, angles and relative directions between the bats
% [DistVec, AngleVec, ~] = FindDistAngle(xBatTx, yBatTx, BatTxDirection, xBatRx, yBatRx, BatRxDirection);
RxInd = find( BATTrantmitter.OtherBatsPolarCoordinates(CurrPulseNum).TargetsID == BatRxNum );
DistVec = BATTrantmitter.OtherBatsPolarCoordinates(CurrPulseNum).Distances(RxInd);
AngleVec = BATTrantmitter.OtherBatsPolarCoordinates(CurrPulseNum).TargetAngle(RxInd);
Lambda = mSoundVelocity ./ abs(TransmitFreqs*1e3); % m
%%
switch CoordinationFlag
    case 'meters'
        xyResolution = 1;
    case 'xy'
        xyResolution= AllParams.SimParams.xyResolution;
end %switch CoordinationFlag

%%
%%%%% Caluculating Direct Interferences %%%%%%%

mDistances = DistVec*xyResolution;
RelativeTransmitAngles = BatTxDirection -  AngleVec; 
% RelativeTransmitAngles = pi*sawtooth(RelativeTransmitAngles-pi);% [-pi,pi]
RelativeTransmitAngles = RelativeTransmitAngles -(2*pi)*floor((RelativeTransmitAngles+pi)/(2*pi));% [-pi,pi]

RelativeRecievedAngles = AngleVec - BatRxDirection + pi;
% RelativeRecievedAngles = pi*sawtooth(RelativeRecievedAngles-pi);% [-pi,pi]
RelativeRecievedAngles = RelativeRecievedAngles -(2*pi)*floor((RelativeRecievedAngles+pi)/(2*pi)); % [-pi,pi]

% try
BatMouthGain = BeamDirectivity(TransmitFreqs, RelativeTransmitAngles, 'Transmit');
% catch
%     pop0 = 'Directivity'
% end

%%%% New 31Oct22
if AllParams.BatSonarParams.DirectivityMode
    BatEarGain = BeamDirectivity(TransmitFreqs, RelativeRecievedAngles, 'Recieved');
else
    BatEarGain = 1;
end

   
%%% THE FORMULA %%%
BatTxAttenuation = AttenuationConst .* BatMouthGain .* BatEarGain .* ...
    10.^(-1.*AlphaAtmAtt./10.*(mDistances-0.1)) .* (Lambda.^2  ./ (4*pi.*mDistances).^2 ); %.*...
        % exp(1i*2*pi.*mDistances./Lambda);
TranmittedPulsePower = 10.^(TranmittedPulsePower/10); % from DB to spl
DirectInterferncePower = TranmittedPulsePower.*BatTxAttenuation;

%%% The times %%%
DelayTimes_precise = mDistances ./ mSoundVelocity ./ SampleTime;
DirectInterfernceTimes_precise = min(Bat1TransmitPulseTimes + DelayTimes_precise, MaxTime);
DelayTimes = round(DelayTimes_precise);   %  round((mDistances ./ mSoundVelocity) ./ SampleTime);
DirectInterfernceTimes = min(Bat1TransmitPulseTimes +DelayTimes, MaxTime) ; % Protecting EndTime

 %% Recontruct the Acoustic Signal in BATReciever for the Direct Interdfrence Vec
    %%%% Nov2021
    if AllParams.SimParams.AcousticsCalcFlag
        FsAcoustic = AllParams.SimParams.FsAcoustic;
        nSamples = numel(BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).TxAcousticSig);
        idx1 = round((TransmitStartTime + DelayTimes_precise) * FsAcoustic * SampleTime);
        EchoDetailed = struct('EchosnTimesFull', DirectInterfernceTimes, 'EchoAttenuationFull', BatTxAttenuation);
        
        Curr_Inter_Acoustic = ReconstractAcousticEcho(EchoDetailed, ...
            BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).TxAcousticSig, SampleTime, nSamples);
        
        BATReciever.AcousticSig_ConspsCalls = Add_currEcho_to_full_vec( ...
            BATReciever.AcousticSig_ConspsCalls , Curr_Inter_Acoustic, idx1, nSamples);
    else
        BATReciever.AcousticSig_ConspsCalls = [];
    end % if AllParams.SimParams.AcousticsCalcFlag
    
%% Removed Nov2021 - Acoustics
%%% FilterBank  %%%%%
% FilterBankFlag = 0;
% FilterBank_jam_struct = struct();
% 
% if strcmp( AllParams.BatSonarParams.ReceiverType, 'FilterBank')
%     FilterBankFlag = 1;
%     fs_ts =  AllParams.BatSonarParams.FilterBank_Fs * AllParams.SimParams.SampleTime;
%     start_us_time = BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).StartPulseTime*fs_ts;
%     %%% remove to AllParams
%     max_delay = 0.03; % sec. equivalamt to 5 meter
%     max_samples = ceil(max_delay*AllParams.BatSonarParams.FilterBank_Fs);
%     max_samples = FilterBank.time_to_analyze_us;
%     
%     max_us_time = min( start_us_time + max_samples, ...
%         MaxTime*fs_ts);
%     
%      % short buzz signals
%         tx_freqs = TransmitFreqs;
%         rx_powers = DirectInterferncePower;
%         if numel(TransmitFreqs) < 2
%             tx_freqs =  BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).ChirpMinMaxFreqs([2,1]); %Hz
%             rx_powers = rx_powers*ones(1,2);
%         end % if numel
%     
%     % init output mat
%     
%     jam_mat = zeros(numel(FilterBank.fft_freqs), max_us_time- start_us_time+1);
%     
%     % update matrix from direct interference
% 
%     % updated_mat
%     [jam_mat, u_fft_idx ] = calc_jam_raw_mat(...
%         jam_mat, DirectInterfernceTimes_precise, tx_freqs, rx_powers, start_us_time,   FilterBank.fft_freqs, fs_ts);
% %     jam_mat = updated_mat; 
%     %   XXXX
% %     jam_mat = jam_mat + updated_mat; 
%     % XXXX
%    
% end % if FilterBank

%%
%%%%% Caluculating Interferences From Other Bats' Echoes ('Phantom Echoes')%%%%%%%

%%% INIT
MaxInterVecSize= AllParams.SimParams.TotalPreysNumber * PulseDuration; 
FullEchosInetfernceTimes = zeros(1,MaxInterVecSize);
FullEchosInetfernceTimes_precise = zeros(1,MaxInterVecSize);
FullEchosInetferncePower = zeros(1,MaxInterVecSize);
FullEchosInetfernceFreq = zeros(1,MaxInterVecSize);
FullEchosInetfernceTxBatNum = zeros(1,MaxInterVecSize);
FullEchosInetfernceTxBatNumMaxRxPowerDB = zeros(1,MaxInterVecSize);
FullEchosInetfernceTxBatNumDistance2PhantomReal = zeros(1,MaxInterVecSize);
FullEchosInetfernceTxBatNumAngle2Phantom = zeros(1,MaxInterVecSize);

NumOfTargets = AllParams.SimParams.TotalPreysNumber;

if  AllParams.BatSonarParams.BiStatMode
    BiStatDetectionSrtct = BATSstrct(BatRxNum).BiStatDetection;
else
    BiStatDetectionSrtct = [];
end % AllParams.BatSonarParams.BiStatMode

if  AllParams.BatSonarParams.PhantomEchoesFromConsFlag 
% %     PhantomEchoesStrct = BATSstrct(BatRxNum).PhantomEchoesStrct;
    PhantomEchoesStrct(NumOfTargets) = struct( ...
        'PreyNum', [], ...
        'TxBatNum',[], ...
        'SignalPhase',[], ... % 'Search', 'Approach', 'Buzz'
        'StartnTime', [], ...
        'EndnTime', [] , ...
        'Freqs', [], ...
        'RxPowers', [], ... 
        'MaxRxPowerDB', [],...
        'Distance2PhantomReal', [],...
        'Angle2Phantom',[] ...
    );
else
    PhantomEchoesStrct = [];
end % AllParams.BatSonarParams.BiStatMode


%%% Setting the parameters
TargetArea = AllParams.BatSonarParams.TargetWingLength.^2*pi;
% TargetRCS = 0.5*TargetArea;

PrevTargetsEchosRxPowers = [];
PrevTargetsEchosRxTimes = [];


%%% The distance  and angle from the transmitter to each target at the
%%% Transmitting times
DistanceTxBat2Target = BATTrantmitter.PreysPolarCoordinate.Distances .* xyResolution;
AngleTxBat2Target = BATTrantmitter.PreysPolarCoordinate.TargetAngle ;
%     DistanceTxBat2Target = BATTrantmitter.Vector2Preys(kTarget).Dist2Prey(Bat1TransmitPulseTimes);
%     AngleTxBat2Target = BATTrantmitter.Vector2Preys(kTarget).Angle2Prey(Bat1TransmitPulseTimes);
%%% The Power
CurrTranmittedPulsePower =  TranmittedPulsePower;

%%% The Dist and Angle of the recived bat from target
%%%(remark: It Calculated at the preivous rxBat Pulse)

% for the first maneuvers , not all bats started to fly
if ~isempty(BATReciever.PreysPolarCoordinate)
    
    DistanceTarget2RxBat =  BATReciever.PreysPolarCoordinate.Distances .* xyResolution;
    AngleTarget2RxBat = BATReciever.PreysPolarCoordinate.TargetAngle;
    
else % % if ~isempty(BATReciever.PreysPolarCoordinate)
    DistanceTarget2RxBat = zeros(1,AllParams.SimParams.TotalPreysNumber);
    AngleTarget2RxBat = zeros(1,AllParams.SimParams.TotalPreysNumber);
    for nPrey = 1:AllParams.SimParams.TotalPreysNumber
        xnPrey = PREY(nPrey).xPrey(TransmitStartTime);
        ynPrey = PREY(nPrey).yPrey(TransmitStartTime);
        nPreyDirection = PREY(nPrey).PreyTeta(TransmitStartTime);
        [DistanceTarget2RxBat(nPrey), AngleTarget2RxBat(nPrey), ~] = ...
            FindDistAngle(xBatRx, yBatRx, BatRxDirection, xnPrey, ynPrey, nPreyDirection);
    end % for nPrey
    
    DistanceTarget2RxBat =DistanceTarget2RxBat * xyResolution;
end % if ~isempty(BATReciever.PreysPolarCoordinate)
for kTarget = 1:NumOfTargets    
    %%% The Dirctivety of Tx and Rx
    TxGain = BeamDirectivity(TransmitFreqs, AngleTxBat2Target(kTarget), 'Transmit');
    RxGain = BeamDirectivity(TransmitFreqs, AngleTarget2RxBat(kTarget), 'Recieved');
    
    %%% Rx Power from Sonar Equation %%%
    SigmaRCS = 1;
%     TargetRCS = nan(size(Lambda));

    % new model for RCS -Nov2022
    TargetRCS = TargetArea*ones(size(Lambda));
    rTarget =  sqrt(TargetArea/(2*pi));
    ixSmall = rTarget./Lambda < 1.5;
%     ixBig = rTarget./Lambda > 4;
    TargetRCS(ixSmall) = SigmaRCS*4*pi*TargetArea^2 ./ Lambda(ixSmall).^2; % New Nov22
%     TargetRCS(~ixSmall & ~ixBig) = SigmaRCS*2*pi*TargetArea; % New Nov22
%     TargetRCS(ixBig) = TargetArea; % New Nov22

    DistTxBat2CurrTarget = DistanceTxBat2Target(kTarget);
    DistCurrTarget2RxBat = DistanceTarget2RxBat(kTarget);
    
    Attenuation = AttenuationConst .* TxGain .* RxGain .* Lambda.^2 .* TargetRCS .* ...
        10.^(-AlphaAtmAtt./10.*(DistTxBat2CurrTarget + DistCurrTarget2RxBat-0.2)) ./ ...
        ( (4*pi).^3 .* DistTxBat2CurrTarget.^2 .* DistCurrTarget2RxBat.^2 ); % .* ...
        % exp(1i*4*pi.*(DistTxBat2CurrTarget + DistCurrTarget2RxBat)./Lambda);
    CurrentEchosInterPower = CurrTranmittedPulsePower .* Attenuation;
    
    %%% Recieving Times
    EchoDelayTimes_precise = (DistTxBat2CurrTarget + DistCurrTarget2RxBat) ./ mSoundVelocity  ./ SampleTime;
    CurrentEchoDelayTimes_precise = min(Bat1TransmitPulseTimes+ EchoDelayTimes_precise, MaxTime);
    EchoDelayTimes = round(EchoDelayTimes_precise);     % round( (DistTxBat2CurrTarget + DistCurrTarget2RxBat) ./ mSoundVelocity  ./ SampleTime);
    CurrentEchosInterTimes = min(Bat1TransmitPulseTimes + EchoDelayTimes , MaxTime);
    
    %%% Concatenated Vectors
    CurrIndex = ((kTarget-1)*PulseDuration+1):(kTarget*PulseDuration);
    FullEchosInetfernceTimes(CurrIndex) = CurrentEchosInterTimes;
    FullEchosInetfernceTimes_precise(CurrIndex) = CurrentEchoDelayTimes_precise;
    FullEchosInetferncePower(CurrIndex) = CurrentEchosInterPower;
    FullEchosInetfernceFreq(CurrIndex) = TransmitFreqs;
    
    %% FilterBank up-smpling
%     %%%% Removed Nov2021
%     if FilterBankFlag
%         % calvulate only phisible echoes  
%         orig_time_vec = FullEchosInetfernceTimes_precise(CurrIndex);
%         if orig_time_vec(1)*fs_ts < max_us_time
%             % short buzz signals
%             rx_powers = CurrentEchosInterPower;
%             if numel(CurrentEchosInterPower) < 2
%                 rx_powers = rx_powers*ones(1,2);
%             end % if numel
%             [jam_mat, ~ ] = calc_jam_raw_mat(...
%                 jam_mat, orig_time_vec, tx_freqs, rx_powers, start_us_time, FilterBank.fft_freqs, fs_ts);
%             
% %             jam_mat = updated_mat;
% 
% %             jam_mat = jam_mat+ updated_mat;
%         end
%         
%         
%     end % if FilterBank
    %% end Filter Bank
    %%
    % Concatenated struct of bistast
    if  AllParams.BatSonarParams.BiStatMode
        Idx = find( CurrentEchosInterPower >=  ...
            BATSstrct(BatRxNum).BiStatDetection.Prey(kTarget).RxPower(CurrentEchosInterTimes));
        IdxUpdate = CurrentEchosInterTimes(Idx);
        BiStatDetectionSrtct.Prey(kTarget).nTimes(IdxUpdate) = IdxUpdate;
        BiStatDetectionSrtct.Prey(kTarget).RxPower(IdxUpdate) = CurrentEchosInterPower(Idx);
        BiStatDetectionSrtct.Prey(kTarget).Freqs(IdxUpdate) = TransmitFreqs(Idx);
        BiStatDetectionSrtct.Prey(kTarget).Dist2Prey(IdxUpdate) = DistCurrTarget2RxBat;
        BiStatDetectionSrtct.Prey(kTarget).Angle2Prey(IdxUpdate) = AngleTarget2RxBat(kTarget);
        BiStatDetectionSrtct.Prey(kTarget).TxBatNum(IdxUpdate) = BATTxNum;
        
        BiStatDetectionSrtct.AllPreyRxPowerMat(kTarget,IdxUpdate) = CurrentEchosInterPower(Idx);
        BiStatDetectionSrtct.AllPFreqMat(kTarget,IdxUpdate) = TransmitFreqs(Idx);
        
    end % if  AllParams.BatSonarParams.BiStatMode
    % update the out put only if current echo-power is greater than
    % pervioues powers in the recption times
    
    %% Output struct for phantom echoes
    if  AllParams.BatSonarParams.PhantomEchoesFromConsFlag
            PrevPulseNum =max(1,CurrPulseNum-1);
            
            PhantomEchoesStrct(kTarget).PreyNum = kTarget ;
            PhantomEchoesStrct(kTarget).TxBatNum = BATTxNum ; 
            PhantomEchoesStrct(kTarget).SignalPhase = BATSstrct(BATTxNum).ManueverCmdStruct(PrevPulseNum).ManueverStage;
            PhantomEchoesStrct(kTarget).StartnTime = CurrentEchosInterTimes(1);
            PhantomEchoesStrct(kTarget).EndnTime = CurrentEchosInterTimes(end);
            PhantomEchoesStrct(kTarget).Freqs = TransmitFreqs;
            PhantomEchoesStrct(kTarget).RxPowersDB = 10*log10(abs(CurrentEchosInterPower));
            PhantomEchoesStrct(kTarget).MaxRxPowerDB = 10*log10(abs(max(CurrentEchosInterPower)));
            PhantomEchoesStrct(kTarget).Distance2PhantomReal = DistCurrTarget2RxBat;
            PhantomEchoesStrct(kTarget).Angle2Phantom = AngleTarget2RxBat(kTarget);
%%        
    end % if  AllParams.BatSonarParams.PhantomEchoesFromConsFlag
    
    %% Recontruct the Acoustic Signal in BATReciever Interdfrence Vec
    %%%% Nov2021
    if AllParams.SimParams.AcousticsCalcFlag
        FsAcoustic = AllParams.SimParams.FsAcoustic;
        nSamples = numel(BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).TxAcousticSig);
        idx1 = round((TransmitStartTime + EchoDelayTimes_precise) * FsAcoustic * SampleTime);
        EchoDetailed = struct('EchosnTimesFull', CurrentEchosInterTimes, 'EchoAttenuationFull', Attenuation);
        
        Curr_Inter_Acoustic = ReconstractAcousticEcho(EchoDetailed, ...
            BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).TxAcousticSig, SampleTime, nSamples);
        
        BATReciever.AcousticSig_ConspsCalls = Add_currEcho_to_full_vec( ...
            BATReciever.AcousticSig_ConspsCalls , Curr_Inter_Acoustic, idx1, nSamples);
    else
        BATReciever.AcousticSig_ConspsCalls= [];
    end % if AllParams.SimParams.AcousticsCalcFlag
end % for kTarget 
    
%%% Checking for times equal between Interferences and sum it %%%
AllInterferncesTimes = [DirectInterfernceTimes, FullEchosInetfernceTimes];
AllInterferncesTimes_precise = [DirectInterfernceTimes_precise, FullEchosInetfernceTimes_precise];
AllInterferncesPower = [DirectInterferncePower, FullEchosInetferncePower];
AllInterferncesFreq = [TransmitFreqs, FullEchosInetfernceFreq];

UniqueTimes = unique(AllInterferncesTimes);
PowerSum = zeros(size(UniqueTimes));
% ALLInterStruct(length(UniqueTimes)) = ...
%     struct('Times',[],'Freqs',[], 'Power', []);

ALLInterStruct = BATReciever.InterferenceFullStruct ;

for kk = 1:length(UniqueTimes)
    CurrTimes = [AllInterferncesTimes == UniqueTimes(kk)];
    PowerSum(kk) = sum( AllInterferncesPower(CurrTimes) );
    try
        ALLInterStruct(UniqueTimes(kk)).Times = UniqueTimes(kk);
        ALLInterStruct(UniqueTimes(kk)).Times_precise = ...
            vertcat(ALLInterStruct(UniqueTimes(kk)).Times_precise, AllInterferncesTimes_precise(CurrTimes)');
        
        %         ALLInterStruct(UniqueTimes(kk)).Freqs = [ALLInterStruct(kk).Freqs, AllInterferncesFreq(CurrTimes)];
        %         ALLInterStruct(UniqueTimes(kk)).Power = [ALLInterStruct(kk).Power, AllInterferncesPower(CurrTimes)];
        
        ALLInterStruct(UniqueTimes(kk)).Freqs =  vertcat(ALLInterStruct(UniqueTimes(kk)).Freqs, AllInterferncesFreq(CurrTimes)');
        ALLInterStruct(UniqueTimes(kk)).Power =  vertcat(ALLInterStruct(UniqueTimes(kk)).Power, AllInterferncesPower(CurrTimes)');
        
        %         ALLInterStruct(UniqueTimes(kk)).BATTxNum = ...
        %             [ALLInterStruct(kk).BATTxNum, BATTxNum*ones(size(AllInterferncesFreq(CurrTimes)))];
        %         ALLInterStruct(UniqueTimes(kk)).BatTxPulseNum = ...
        %             [ALLInterStruct(kk).BatTxPulseNum, CurrPulseNum*ones(size(AllInterferncesFreq(CurrTimes)))];
        
        ALLInterStruct(UniqueTimes(kk)).BATTxNum = ...
             vertcat(ALLInterStruct(UniqueTimes(kk)).BATTxNum, BATTxNum*ones(size(AllInterferncesFreq(CurrTimes)))');
        ALLInterStruct(UniqueTimes(kk)).BatTxPulseNum = ...
             vertcat(ALLInterStruct(UniqueTimes(kk)).BatTxPulseNum, CurrPulseNum*ones(size(AllInterferncesFreq(CurrTimes)))');
    catch
        hopa ='strct inter'

    end% try
end % for kk

%% Output
AllInetfernceTimes  = UniqueTimes;
AllInetferncePower  = PowerSum;

% update the Acoustic InterferenceSignal
AcousticSig_ConspsCalls = BATReciever.AcousticSig_ConspsCalls;

% %%% FilterBank OutPut - Removed Nov2021
% if FilterBankFlag
%     FilterBank_jam_struct = struct( ...
%         'jam_raw_mat', jam_mat, ...
%         'unique_fft_idx', u_fft_idx, ...
%         'start_update_time', start_us_time, ...
%         'end_update_time', max_us_time...
%         );
% end % if FilterBankFllag

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Power] = AddPowerIfTimeExists(TimesVec,PowerVec,CheckTime)
    Ind = find(TimesVec == CheckTime);
    if ~isempty(Ind)
        Power = PowerVec(Ind);
    else
        Power = 0;
    end
    %%% dealng with the end of the time vector (mulpilie echos at nTime=
    %%% end)
    if length(Power) > 1 
        Power = Power(1,1);
    end % if length
end