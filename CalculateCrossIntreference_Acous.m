function [AcousticSig_Consps, TransmitFreqs, rxTable] = ...
    CalculateCrossIntreference_Acous...
    (BATSstrct, BATTxNum, BatRxNum, PREY, CoordinationFlag, AllParams, CurrPulseNum, SigType, EarSide)
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
%
%           SigType = 'ConspsCalls' | 'ConspsEchoes' 
%
% OutPuts:  
%           Bat2RecievedPulseTimes - the nTimes vector of the recieved interference
%           Bat2Attenuation - the vector of the received powers 
%           EchosIntRecievedPower, EchosIntRecievedTimes - the power and
%               times of the echos of targets from the interferer bat
%           Freqs - vector of the frequncy recieved
%           rxTabe = a table of the received signals oove TH level, the
%           table consists: time, callNum, distance, direction ofArraival,
%           receivedLevel, TxBat, and echoType
%

% 15/04/19 updated by omer : Adding bi-stat detection output

%%% 7/11/21 
%%% The Building of the Acoustic Signal of the Interference is
%%% executerd for each kRxBat in the function CalculateCrossIntreferenceOnLine
%%%%% NOV2021

%%% 11Sep22
% Acoustiocs Double reflectiones

%%%%% 24/1/23
% Add ear Side to Directivity
if nargin < 9
   EarSide = 'None' ;
end

%%
%%%% Inputs

BATTrantmitter = BATSstrct(BATTxNum);
BATReciever    =  BATSstrct(BatRxNum);
%%% CONSTS and Parameters from Inputs
mSoundVelocity = AllParams.SimParams.SoundV0; %m/sec
SampleTime = AllParams.SimParams.SampleTime;
xyResolution = AllParams.SimParams.xyResolution;
MaxTime = AllParams.SimParams.SimulationTime ./ SampleTime ;
% TranmittedPulsePower = 10^(BATTrantmitter.PulsePower/10);
AttenuationConst = 1; % Just parameter

switch EarSide
    case 'None'
        relRxAngle = 0;
    case 'Left'
        relRxAngle = AllParams.BatSonarParams.EarAngle;
    case 'Right'
        relRxAngle = -AllParams.BatSonarParams.EarAngle;
end % switch
% Nov2021

%%   
%%%% init
rxTable = struct();

%%% the tranmitted time, power and freqs
TransmitStartTime = BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).StartPulseTime;
Bat1TransmitPulseTimes = TransmitStartTime:(TransmitStartTime+ BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).PulseDuration-1); 
TransmitFreqs = BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).PulseFreqsCommands;
TranmittedPulsePower = BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).PulsePower;
TranmittedPulsePowerLin = 10.^(TranmittedPulsePower /10);
PulseDuration = length(TransmitFreqs);
%%
AlphaAtmAtt = 0.038*TransmitFreqs - 0.3; %  from table, @20Cel.50% humidity Atmosphere attenuation[m^(-1)] ...
                    %  1dB loose/meter % 0.3 for 40kHa (at 30C,70% humidty],alpha =  0.1 for 28kH
Lambda = mSoundVelocity ./ abs(TransmitFreqs*1e3); % m

%%
%%% the position of the bats when Bat1 transmits
xBatTx = deal(BATTrantmitter.xBati(TransmitStartTime));
yBatTx = deal(BATTrantmitter.yBati(TransmitStartTime));
BatTxDirection = deal(BATTrantmitter.Teta(TransmitStartTime));
xBatRx = deal(BATReciever.xBati(TransmitStartTime));
yBatRx = deal(BATReciever.yBati(TransmitStartTime));
BatRxDirection = deal(BATReciever.Teta(TransmitStartTime));


%%
switch CoordinationFlag
    case 'meters'
        xyResolution = 1;
    case 'xy'
        xyResolution= AllParams.SimParams.xyResolution;
end %switch CoordinationFlag

switch SigType
    case 'ConspsCalls'
        %% Caluculating Direct Interferences %%%%%%%
        
        % The Full acoustic signal:
        AcousticSig_Consps = BATReciever.AcousticSig_ConspsCalls;
        
        % Calculate RxLvl by Cross Sonar equation
        CalcCrossSonarFlg = 0;
        %%% the distances, angles and relative directions between the bats
        % [DistVec, AngleVec, ~] = FindDistAngle(xBatTx, yBatTx, BatTxDirection, xBatRx, yBatRx, BatRxDirection);
        RxInd = find( BATTrantmitter.OtherBatsPolarCoordinates(CurrPulseNum).TargetsID == BatRxNum );
        DistVec = BATTrantmitter.OtherBatsPolarCoordinates(CurrPulseNum).Distances(RxInd);
        AngleVec = BATTrantmitter.OtherBatsPolarCoordinates(CurrPulseNum).TargetAngle(RxInd);

        mDistances = DistVec*xyResolution;
        RelativeTransmitAngles = BatTxDirection -  AngleVec;
        % RelativeTransmitAngles = pi*sawtooth(RelativeTransmitAngles-pi);% [-pi,pi]
        RelativeTransmitAngles = RelativeTransmitAngles -(2*pi)*floor((RelativeTransmitAngles+pi)/(2*pi));% [-pi,pi]

        RelativeRecievedAngles = AngleVec - (BatRxDirection+relRxAngle)  + pi;
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
        InterferncesPower = TranmittedPulsePower.*BatTxAttenuation;

        %%% The times %%%
        DelayTimes_precise = mDistances ./ mSoundVelocity ./ SampleTime;
%         InterfernceTimes_precise = min(Bat1TransmitPulseTimes + DelayTimes_precise, MaxTime);
        DelayTimes = round(DelayTimes_precise);   %  round((mDistances ./ mSoundVelocity) ./ SampleTime);
        rxTimes = min(Bat1TransmitPulseTimes +DelayTimes, MaxTime) ; % Protecting EndTime
        
        %% Recontruct the Acoustic Signal in BATReciever for the Direct Interdfrence Vec
        %%%% Nov2021
        if AllParams.SimParams.AcousticsCalcFlag
            FsAcoustic = AllParams.SimParams.FsAcoustic;
            nSamples = numel(BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).TxAcousticSig);
            idx1 = round((TransmitStartTime + DelayTimes_precise) * FsAcoustic * SampleTime);
            EchoDetailed = struct('EchosnTimesFull', rxTimes, 'EchoAttenuationFull', BatTxAttenuation);

            Curr_Inter_Acoustic = ReconstractAcousticEcho(EchoDetailed, ...
                BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).TxAcousticSig, SampleTime, nSamples);

%             AcousticSig_Consps = Add_currEcho_to_full_vec( ...
%                 AcousticSig_Consps , Curr_Inter_Acoustic, idx1, nSamples);

            %% %%%% XXX YYY XXX
            % the indices for the assignment
            Acoustic_idx = idx1 : idx1 + nSamples-1;
             AcousticNumel = numel(AcousticSig_Consps);

            % check for the end of simulation
            if idx1 < AcousticNumel

                if Acoustic_idx(end) > AcousticNumel
                    Acoustic_idx = Acoustic_idx(1):AcousticNumel;
                    Curr_Inter_Acoustic = Curr_Inter_Acoustic(1:numel(Acoustic_idx));
                end % if Acoutic_idx(end) > AcousticNumel
                %%% add the result to the acoustic vector
                %     full_vec(Acoustic_idx) = prev_full(Acoustic_idx) + curr_vec;
                AcousticSig_Consps(Acoustic_idx) = AcousticSig_Consps(Acoustic_idx) + Curr_Inter_Acoustic;

            else % if idx1 < AcousticNumel
                % do nothing
            end % if idx1 < AcousticNumel

            %% %%%% XXX YYY XXX
        else
            AcousticSig_Consps = [];
        end % if AllParams.SimParams.AcousticsCalcFlag

        %%   % times, Distances and angles for rxTable- New June2024
        try
            N = numel(mDistances);
            rxTable.time               = TransmitStartTime*SampleTime + mDistances(:)./mSoundVelocity;
            rxTable.CallNum            = CurrPulseNum*ones(N,1);
            rxTable.distance           = mDistances(:);
            rxTable.directionOfArrival = RelativeRecievedAngles(:);
            rxTable.receivedLevel      = max(InterferncesPower(:));
            rxTable.TxBat              = BATTxNum*ones(N,1);
            rxTable.echoType           = repmat(categorical("consCall"), N,1);
            % taregt index is unique value for each rx signal
            rxTable.targetIndex        = randi(1e12,[N,1]); % rxTable.TxBat + linspace(0.1, 0.9 ,N);  % taregt index is unique value for each rx signal 
            rxTable = struct2table(rxTable);
        catch
            warning('Calc rxTable Probllem')
        end
        
    case 'ConspsPreyEchoes'
        %% Caluculating Interferences From Other Bats' Echoes ('Phantom Echoes')%%%%%%%
        
        % The Full acoustic signal:
        AcousticSig_Consps = BATReciever.AcousticSig_ConspsPreyEchoes;
        
        % Calculate RxLvl by Cross Sonar equation
        CalcCrossSonarFlg = 1;
        %%% INIT
%         MaxInterVecSize= AllParams.SimParams.TotalPreysNumber * PulseDuration;
%         InterfernceTimes = zeros(1,MaxInterVecSize);
%         InterfernceTimes_precise = zeros(1,MaxInterVecSize);
%         InterferncesPower = zeros(1,MaxInterVecSize);
%         FullEchosInetfernceFreq = zeros(1,MaxInterVecSize);
        
        NumOfTargets = AllParams.SimParams.TotalPreysNumber;

        %%% Setting the parameters
        TargetArea = AllParams.BatSonarParams.TargetWingLength.^2*pi;
        % TargetRCS = 0.5*TargetArea;

        %%% The distance  and angle from the transmitter to each target at the
        %%% Transmitting times
        DistanceTxBat2Target = BATTrantmitter.PreysPolarCoordinate.Distances .* xyResolution;
        AngleTxBat2Target    = wrapToPi(BATTrantmitter.PreysPolarCoordinate.TargetAngle);
        %     DistanceTxBat2Target = BATTrantmitter.Vector2Preys(kTarget).Dist2Prey(Bat1TransmitPulseTimes);
        %     AngleTxBat2Target = BATTrantmitter.Vector2Preys(kTarget).Angle2Prey(Bat1TransmitPulseTimes);
        %%% The Power
        CurrTranmittedPulsePower =  TranmittedPulsePower;

        %%% The Dist and Angle of the recived bat from target
        %%%(remark: It is Calculated at the preivous rxBat Pulse)

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

    case {'ConspsObsEchoes', 'LeftObs', 'RightObs'}
     %% Caluculating Interferences From Other Bats' Echoes ('Phantom Echoes')%%%%%%%
        
        % The Full acoustic signal:
        fieldName = strcat('AcousticSig_', SigType);
        AcousticSig_Consps = BATReciever.(fieldName);
        
        % Calculate RxLvl by Cross Sonar equation
        CalcCrossSonarFlg = 1; 

        % The relevant clutter is in the struct BATTrantmitter.ObsInBeamStruct
        ObsInBeamStruct = BATTrantmitter.ObsInBeamStruct(CurrPulseNum);

        NumOfTargets = numel(ObsInBeamStruct.Distances);

        %%% Setting the parameters
        TargetArea = AllParams.BatSonarParams.ClutterPointReflecetingArea .* 10.^(AllParams.BatSonarParams.ClutterGainDB/10);

        %%% The distance  and angle from the transmitter to each target at the
        %%% Transmitting times
        DistanceTxBat2Target = ObsInBeamStruct.Distances .* xyResolution;
        AngleTxBat2Target    = wrapToPi(ObsInBeamStruct.TargetAngle);

        xObs = ObsInBeamStruct.xOBSTC;
        yObs = ObsInBeamStruct.yOBSTC;
        
        %  calculate the polar coordinations relative to Rx Bat
        dY = yObs - yBatRx  ;
        dX = xObs - xBatRx ;
        DistanceTarget2RxBat = sqrt(dY.^2 +  dX.^2);
        DistanceTarget2RxBat = DistanceTarget2RxBat * xyResolution;
        AngleTarget2RxBat    = wrapToPi(atan2(dY, dX) - (BatRxDirection+relRxAngle));
        
        %%% The Power
        CurrTranmittedPulsePower =  TranmittedPulsePower;

    case 'ConspsConspsEchoes' 
        %% Caluculating Interferences From Other Bats' Echoes ('Phantom Echoes')%%%%%%%

        % The Full acoustic signal:
        AcousticSig_Consps = BATReciever.AcousticSig_ConspsConspsEchoes;
        % Calculate RxLvl by Cross Sonar equation
        CalcCrossSonarFlg = 1;
        % The relevant consps is in the struct BATTrantmitter.ObsInBeamStruct
        ConspsInBeamStruct = BATTrantmitter.OtherBatsPolarCoordinates(CurrPulseNum);
        
        % Remove own reciever from conspecifics
        ixRx = ConspsInBeamStruct.TargetsID == BatRxNum;
        NumOfTargets = numel(ConspsInBeamStruct.Distances);
        for s = fieldnames(ConspsInBeamStruct)'
            if numel(ConspsInBeamStruct.(s{:})) == NumOfTargets
              ConspsInBeamStruct.(s{:})(ixRx)  = [];
            end
        end % for s
        NumOfTargets = numel(ConspsInBeamStruct.Distances);
        ConspsInBeamStruct.NumOfTargets = NumOfTargets;
        %%% Setting the parameters
        TargetArea = AllParams.BatSonarParams.BatWingLength^2 * pi;
       

        %%% The distance  and angle from the transmitter to each target at the
        %%% Transmitting times
        DistanceTxBat2Target = ConspsInBeamStruct.Distances .* xyResolution;
        AngleTxBat2Target    = wrapToPi(ConspsInBeamStruct.TargetAngle);
        %%% The distance  and angle from the receiver to each target
        DistanceTarget2RxBat = nan(size(DistanceTxBat2Target));
        AngleTarget2RxBat    = nan(size(DistanceTxBat2Target));
end % switch


%% Calculate the Cross Sonar equation for Consps Echoes

SigmaRCS = 1;
if CalcCrossSonarFlg
    allAtten = nan(NumOfTargets,1);

    for kTarget = 1:NumOfTargets
        %%% find the position of the conspecifics in the current pulse start
        if strcmp('ConspsConspsEchoes', SigType )
            targetIx = ConspsInBeamStruct.TargetsID(kTarget); 
            xC = BATSstrct(targetIx).xBati(TransmitStartTime);
            yC = BATSstrct(targetIx).yBati(TransmitStartTime);
            %  calculate the polar coordinations relative to Rx Bat
            dY = yC - yBatRx  ;
            dX = xC - xBatRx ;
            DistanceTarget2RxBat(kTarget) = sqrt(dY.^2 +  dX.^2) * xyResolution;
            AngleTarget2RxBat(kTarget)    = wrapToPi(atan2(dY, dX) - (BatRxDirection+relRxAngle)) ; % Add  Left/ Right ear relative direction
        end % if strcmp('ConspsConspsEchoes', SigType )
        %%% The Dirctivety of Tx and Rx
        TxGain = BeamDirectivity(TransmitFreqs, AngleTxBat2Target(kTarget), 'Transmit');
        %%%% Add  Left/ Right ear relative direction
%         AngleTarget2RxBat(kTarget) = AngleTarget2RxBat(kTarget) + relRxAngle; 
        %%%    
        RxGain = BeamDirectivity(TransmitFreqs, AngleTarget2RxBat(kTarget), 'Recieved');

        %%% Rx Power from Sonar Equation %%% Fixed 30Nov22
        % TargetRCS = 4*pi*TargetArea^2./Lambda.^2;
        rTarget =  sqrt(TargetArea/pi);
        TargetRCS = TargetArea*ones(size(Lambda));
        ixSmall = rTarget./Lambda < 1.5;
        %     ixBig = rTarget./Lambda > 4;
        TargetRCS(ixSmall) = SigmaRCS*4*pi*TargetArea^2 ./ Lambda(ixSmall).^2; % New Nov22
        
        DistTxBat2CurrTarget = DistanceTxBat2Target(kTarget);
        DistCurrTarget2RxBat = DistanceTarget2RxBat((kTarget));

        % Test Only Target in range (for time performance)
        if (DistTxBat2CurrTarget + DistCurrTarget2RxBat) < AllParams.BatSonarParams.BatDetectionRange
            Attenuation = AttenuationConst .* TxGain .* RxGain .* Lambda.^2 .* TargetRCS .* ...
                10.^(-AlphaAtmAtt./10.*(DistTxBat2CurrTarget + DistCurrTarget2RxBat-0.2)) ./ ...
                ( (4*pi).^3 .* DistTxBat2CurrTarget.^2 .* DistCurrTarget2RxBat.^2 ); % .* ...
            % exp(1i*4*pi.*(DistTxBat2CurrTarget + DistCurrTarget2RxBat)./Lambda);

            %         CurrentEchosInterPower = CurrTranmittedPulsePower .* Attenuation;
            EchoDelayTimes_precise = (DistTxBat2CurrTarget + DistCurrTarget2RxBat) ./ mSoundVelocity  ./ SampleTime;
            EchoDelayTimes = round(EchoDelayTimes_precise);     % round( (DistTxBat2CurrTarget + DistCurrTarget2RxBat) ./ mSoundVelocity  ./ SampleTime);
            CurrentEchosInterTimes = min(Bat1TransmitPulseTimes + EchoDelayTimes , MaxTime);
           allAtten(kTarget) = max(Attenuation);

            %% Recontruct the Acoustic Signal in BATReciever Interdfrence Vec
            %%%% Nov2021
            if AllParams.SimParams.AcousticsCalcFlag
                FsAcoustic = AllParams.SimParams.FsAcoustic;
                nSamples = numel(BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).TxAcousticSig);
                idx1 = round((TransmitStartTime + EchoDelayTimes_precise) * FsAcoustic * SampleTime);
                EchoDetailed = struct('EchosnTimesFull', CurrentEchosInterTimes, 'EchoAttenuationFull', Attenuation);

                Curr_Inter_Acoustic = ReconstractAcousticEcho(EchoDetailed, ...
                    BATTrantmitter.TransmittedPulsesStruct(CurrPulseNum).TxAcousticSig, SampleTime, nSamples);

%                 AcousticSig_Consps = Add_currEcho_to_full_vec( ...
%                     AcousticSig_Consps , Curr_Inter_Acoustic, idx1, nSamples);
                
                %% %%%% XXX YYY XXX
                % the indices for the assignment
                Acoustic_idx = idx1 : idx1 + nSamples-1;
                AcousticNumel = numel(AcousticSig_Consps);

                % check for the end of simulation
                if idx1 < AcousticNumel

                    if Acoustic_idx(end) > AcousticNumel
                        Acoustic_idx = Acoustic_idx(1):AcousticNumel;
                        Curr_Inter_Acoustic = Curr_Inter_Acoustic(1:numel(Acoustic_idx));
                    end % if Acoutic_idx(end) > AcousticNumel
                    %%% add the result to the acoustic vector
                    %     full_vec(Acoustic_idx) = prev_full(Acoustic_idx) + curr_vec;
                    AcousticSig_Consps(Acoustic_idx) = AcousticSig_Consps(Acoustic_idx) + Curr_Inter_Acoustic;

                else % if idx1 < AcousticNumel
                    % do nothing
                end % if idx1 < AcousticNumel

                %% %%%% XXX YYY XXX

            else
                AcousticSig_Consps= [];
            end % if AllParams.SimParams.AcousticsCalcFlag
        
        end % if (DistTxBat2CurrTarget + DistCurrTarget2RxBat) < AllParams.BatSonarParams.BatDetectionRange
    end % for kTarget
end % if CalcCrossSonarFalg

%%   % times, Distances and angles for rxTable- New June2024
% update the rxTable only if it ewas not calculated before
% (ConeseocicicCalls)
try
    if ~istable(rxTable) % && numel(DistanceTarget2RxBat) > 0
        N = numel(DistanceTarget2RxBat);
        rxTable.time               = TransmitStartTime*SampleTime + (DistanceTarget2RxBat(:) + DistanceTxBat2Target(:)) ./ mSoundVelocity ;
        rxTable.CallNum            = CurrPulseNum*ones(N,1);
        rxTable.distance           = DistanceTarget2RxBat(:);
        rxTable.directionOfArrival = AngleTarget2RxBat(:);
        rxTable.receivedLevel      = TranmittedPulsePowerLin * allAtten(:);
        rxTable.TxBat              = BATTxNum*ones(size(rxTable.directionOfArrival));
        rxTable.echoType           = repmat(categorical(string(SigType)), N,1);
        % taregt index is unique value for each rx signal
        rxTable.targetIndex        =  randi(1e12,[N,1]);  % (rxTable.TxBat + linspace(0.1, 0.9 ,N)';  

        rxTable = struct2table(rxTable);
    end
catch
    warning('Calc rxTable Probllem')
end % try 


%%%%% Deleted Sep2022   
%% Checking for times equal between Interferences and sum it %%%


%% Output
% AllInetfernceTimes  = UniqueTimes;
% AllInetferncePower  = PowerSum;

% update the Acoustic InterferenceSignal
% AcousticSig_Consps = BATReciever.AcousticSig_ConspsCalls;


end % function
