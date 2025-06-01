function[BatTxRecievedPulseTimes, BatTxRecievedPulsePower, AllTargetEchosRecievedTimes, AllTargetEchosRecievedTPower, Freqs] = ...
    CalculateCrossIntreference(BATTrantmitter, BATReciever, PREY, CoordinationFlag, AllParams, varargin)

% this function calcultes and returns the interfence that BATTrantmitter causes to
% BATReciever
% using the equation : Pr= Pt*Dt*Dr*(lambda/(4*pi*R))^2 * exp(alpha*R)
% where"    Pr= Pweor recieved
%           Pt= tramsmitted Pwer
%           Dt- Directivity of transmitter (antenna). , function pf angle
%           Dr - Directivity of reciever (antenna), function pf angle
%           lambda = wavelength = SoundVelocity/ Frequency
%           alpha- Atmosphere absorbtion parameter
%           R - distance
% OutPuts:  
%           Bat2RecievedPulseTimes - the nTimes vector of the recieved interference
%           Bat2Attenuation - the vector of the reveived powers 
%           EchosIntRecievedPower, EchosIntRecievedTimes - the power and
%               times of the echos of targets from the interferer bat
%           Freqs - vector of the frequncy recieved
%           I
%

%%%% Inputs

%%% CONSTS and Parameters from Inputs
mSoundVelocity = AllParams.SimParams.SoundV0; %m/sec
SampleTime = AllParams.SimParams.SampleTime;
xyResolution = AllParams.SimParams.xyResolution;
MaxTime = AllParams.SimParams.SimulationTime ./ SampleTime ;
% TranmittedPulsePower = 10^(BATTrantmitter.PulsePower/10);
AttenuationConst = 1; % Just parameter


                    
%%% the tranmitted freqs' zero means no transmittion                    
TransmitFreqs= BATTrantmitter.PulseCurrFreqVec;
[R, Bat1TransmitPulseTimes, Freqs] = find(TransmitFreqs);
	 % Protecting the end time of the vector 

% TransmitPower = BATTrantmitter.PulsePowerVec;
[R, Bat1TransmitPulseTimesPower, TranmittedPulsePower] = find( BATTrantmitter.PulsePowerVec);

TotalNumOfPulses = length(Bat1TransmitPulseTimes);
AlphaAtmAtt = 0.034*Freqs - 0.12; %  from table, @20Cel.50% humidity Atmosphere attenuation[m^(-1)] ...
                    %  1dB loose/meter % 0.3 for 40kHa (at 30C,70% humidty],alpha =  0.1 for 28kH

%%% the position of the bats when Bat1 transmits
xBatTx = BATTrantmitter.xBatPos(Bat1TransmitPulseTimes);
yBatTx = BATTrantmitter.yBatPos(Bat1TransmitPulseTimes);
BatTxDirection = BATTrantmitter.Teta(Bat1TransmitPulseTimes);
xBatRx = BATReciever.xBatPos(Bat1TransmitPulseTimes);
yBatRx = BATReciever.yBatPos(Bat1TransmitPulseTimes);
BatRxDirection = BATReciever.Teta(Bat1TransmitPulseTimes);

%%% the distances, angles and relative directions between the bats
[DistVec, AngleVec, RelativeDirection] = FindDistAngle(xBatTx, yBatTx, BatTxDirection, xBatRx, yBatRx, BatRxDirection);
                                     
Lambda = mSoundVelocity ./ abs(Freqs*1e3); % m

switch CoordinationFlag
    case 'meters'
        mDistances = DistVec;
    case 'xy'
        xyResolution= varargin{1};
        mDistances = DistVec*xyResolution;
end %switch CoordinationFlag

%%%%% Caluculating Direct Interferences %%%%%%%

RelativeTransmitAngles = BatTxDirection -  AngleVec; 
RelativeTransmitAngles = pi*sawtooth(RelativeTransmitAngles-pi);% [-pi,pi]
RelativeRecievedAngles = AngleVec - BatRxDirection + pi;
RelativeRecievedAngles = pi*sawtooth(RelativeRecievedAngles-pi);% [-pi,pi]
try
    BatMouthGain = BeamDirectivity(Freqs, RelativeTransmitAngles, 'Transmit');
catch
    pop0 = 'Directivity'
end
BatEarGain = BeamDirectivity(Freqs, RelativeRecievedAngles, 'Recieved');
   
%%% THE FORMULA %%%
BatTxAttenuation = AttenuationConst .* BatMouthGain .* BatEarGain .* ...
    exp(-1*AlphaAtmAtt.*(mDistances-0.1)) ./ (4*pi.*mDistances).^2 .*...
        exp(1i*2*pi.*mDistances./Lambda);
TranmittedPulsePower = 10.^(TranmittedPulsePower/10); % from DB to spl
BatTxRecievedPulsePower = TranmittedPulsePower.*BatTxAttenuation;

%%% The times %%%
DelayTimes = round((DistVec ./ mSoundVelocity) ./ SampleTime);
BatTxRecievedPulseTimes = min(Bat1TransmitPulseTimes +DelayTimes, MaxTime) ; % Protecting EndTime


%%%%% Caluculating Interferences From Other Bats' Echos%%%%%%%

TargetArea = 0.02^2*2*pi;
TargetRCS = TargetArea;
SizeP = size(PREY);
NumOfTargets = SizeP(1,2);
PrevTargetsEchosRecievedPowers = [];
PrevTargetsEchosRecievedTimes = [];

for kTarget = 1:NumOfTargets
    
    %%% The distance  and angle from the transmitter to each target at the
    %%% Transmitting times
    DistanceTxBat2Target = BATTrantmitter.Vector2Preys(kTarget).Dist2Prey(Bat1TransmitPulseTimes);
    AngleTxBat2Target = BATTrantmitter.Vector2Preys(kTarget).Angle2Prey(Bat1TransmitPulseTimes);
    %%% The Power
    CurrTranmittedPulsePower =  10.^(BATTrantmitter.PulsePowerVec(Bat1TransmitPulseTimes)./10);
    
    %%% The Dist and Angle of the recived bat from target
    DistanceTarget2RxBat =  BATReciever.Vector2Preys(kTarget).Dist2Prey(Bat1TransmitPulseTimes);
    AngleTarget2RxBat =BATReciever.Vector2Preys(kTarget).Angle2Prey(Bat1TransmitPulseTimes);
    
    %%% The Dirctivety of Tx and Rx
    TxGain = BeamDirectivity(Freqs, AngleTxBat2Target, 'Transmit');
    RxGain = BeamDirectivity(Freqs, AngleTarget2RxBat, 'Recieved');
    
    %%% Rx Power from Sonar Equation %%%
    Attenuation = AttenuationConst .* TxGain .* RxGain .* Lambda.^2 .* TargetRCS .* ...
        10.^(-AlphaAtmAtt./10.*(DistanceTxBat2Target + DistanceTarget2RxBat-0.2)) ./ ...
        ( (4*pi).^3 .* DistanceTxBat2Target.^2 .* DistanceTarget2RxBat.^2 ) .* ...
        exp(1i*4*pi.*(DistanceTxBat2Target + DistanceTarget2RxBat)./Lambda);
    CurrentTargetEchosRecievedPower = CurrTranmittedPulsePower .* Attenuation;
    %
    %%% Recieving Times
    EchoDelayTimes = round( (DistanceTxBat2Target + DistanceTarget2RxBat) ./ mSoundVelocity  ./ SampleTime);
    CurrentTargetEchosRecievedTimes = min(Bat1TransmitPulseTimes + EchoDelayTimes , MaxTime);
    
    %%% Checking for times equal between targets and adding power %%%
    AllTargetEchosRecievedTimes = unique([ PrevTargetsEchosRecievedTimes, CurrentTargetEchosRecievedTimes ]);
    NumOfAllTimes = length(AllTargetEchosRecievedTimes);
    % Init Sums
    CurrentTargetPower = zeros(1,NumOfAllTimes);
    PrevTargetPower = zeros(1,NumOfAllTimes);
    AllTargetEchosRecievedTPower = zeros(1,NumOfAllTimes);
    for kEcho = 1: NumOfAllTimes
        TimeToCheck =  AllTargetEchosRecievedTimes(kEcho);
        try
            CurrentTargetPower(kEcho) = AddPowerIfTimeExists( CurrentTargetEchosRecievedTimes, CurrentTargetEchosRecievedPower, TimeToCheck);
            PrevTargetPower(kEcho) = AddPowerIfTimeExists( PrevTargetsEchosRecievedTimes, PrevTargetsEchosRecievedPowers, TimeToCheck);
        catch
            pop ='CurrentTargetEchosRecievedTimes'
        end % try
        AllTargetEchosRecievedTPower(kEcho) = CurrentTargetPower(kEcho)+PrevTargetPower(kEcho);
    end % for kEchos = 1: NumOfAllTimes
    PrevTargetsEchosRecievedTimes = AllTargetEchosRecievedTimes;
    PrevTargetsEchosRecievedPowers = AllTargetEchosRecievedTPower;
end %for kTarget = 1:NumOfTargets

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