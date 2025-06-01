function RxStruct = calcSigRxlvl...
    (FindsStruct, BatOrPreyFlag , CurrentFreqs, nTime, CurrentPulseNum, varargin)
%%
%% NEW for Jamming Moths
% this function calcultes and returns the Rx level the the prey items
% receive the bat's calls or the level recieved by the batrs from the
% Singla Power

% using the equation : Pr= Pt*Dt*Dr*(lambda/(4*pi*R))^2 * exp(alpha*R)
% where:    
%           Pr= Power recieved
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

% 16/04/22 updated by omer : 
%
NumberOfInputs= nargin; 
if NumberOfInputs >= 5
    AllParams =  varargin{1};
    xyResolution = AllParams.SimParams.xyResolution;
    mSoundVelocity = AllParams.SimParams.SoundV0; % Sound Velocity in meters
    SampleTime = AllParams.SimParams.SampleTime;
    NumberofBats = AllParams.SimParams.TotalBatsNumber;
    NumberofPreys = AllParams.SimParams.TotalPreysNumber;
else % if NumberOfInputs >= 5;
    xyResolution = 0.01;
    mSoundVelocity = 343;
    SampleTime = 0.001;
end %if NumberOfInputs >= 5;
if NumberOfInputs == 7
    NumOfEchos =  varargin{2};
    DetectedPreys = varargin{3};
end %if NumberOfInputs == 7

Distances  = FindsStruct.Distances;
Angles     = FindsStruct.TargetAngle;
NumOfEchos = FindsStruct.NumOfTargets;

%%% calculations
    AttenuationConst = 1;
    Lambda = mSoundVelocity ./ abs(CurrentFreqs*1e3); % m
    AlphaAtmAtt = 0.038*CurrentFreqs - 0.3; %  from table, @20Cel.50% humidity Atmosphere attenuation[m^(-1)] ...
                    %  1dB loose/meter % 0.3 for 40kHa (at 30C,70% humidty],alpha =  0.1 for 28kH

%     NumOfEchos = length(FindsStruct.Distance); % How many differnt angles are ecohing...
    RxStruct.NumOfEchos = NumOfEchos;
    RxStruct.PulseDuration = length(CurrentFreqs);
    [SortedDistances, IndexSorting] = sort(Distances);
    RxStruct.TargetIndex = IndexSorting;

%     EchosStruct.EchosTimes = SortedEcosDistances*2./xySoundVelocity ;   % the delays of the echos' sorted from nearest to farset
    mDistances = SortedDistances*xyResolution; % Distnatnces in meters
    Angles = Angles(IndexSorting);
    RxStruct.EchosAngles = Angles ;
    RxStruct.TransmittedPulseTime = nTime;
    RxStruct.EchosTimes = mDistances./mSoundVelocity ./SampleTime + RxStruct.TransmittedPulseTime;
    RxStruct.EchosmDistances = mDistances;
    RxStruct.EchosCenterFreq = median(CurrentFreqs); % the median frequncy for flaculation
    RxStruct.EchosAttenuation = nan(size(SortedDistances));
    RxStruct.TransmittedPulseNum = CurrentPulseNum;
    
    RxStruct.EchoDetailed(max(NumOfEchos,1)) = struct('TargetIndex',[], 'EchoAttenuationFull',[], 'EchosnTimesFull',[], 'EchosFreqsFull',[]);
    for kk = 1:NumOfEchos
        RoundTime = round([RxStruct.EchosTimes(kk)]);
        RxStruct.EchoDetailed(kk).TargetIndex = IndexSorting(kk);
        switch BatOrPreyFlag
            case 'Prey' 
                TxGain = BeamDirectivity(CurrentFreqs, Angles(kk), 'Transmit');
                RxGain = 1; % The moth hearing is omidirectional
            case 'Bat'
                TxGain = 1; % The moth hearing is omidirectional
                RxGain = BeamDirectivity(CurrentFreqs, Angles(kk), 'Recieved');
        end % switch    
        %%% THE FORMULA %%%
        RxStruct.EchoDetailed(kk).EchoAttenuationFull = AttenuationConst .* TxGain .* RxGain .* ...
            10.^(-1.*AlphaAtmAtt./10.*(mDistances(kk)-0.1)) .* (Lambda.^2  ./ (4*pi.*mDistances(kk)).^2 ); %.*...
                % exp(1i*2*pi.*mDistances./Lambda);      
        RxStruct.EchoDetailed(kk).EchosnTimesFull = ...
            RxStruct.TransmittedPulseTime + [RoundTime : ( RoundTime +RxStruct.PulseDuration-1 )];
        RxStruct.EchoDetailed(kk).EchosFreqsFull = CurrentFreqs;
        RxStruct.EchosAttenuation(kk) = max(abs(RxStruct.EchoDetailed(kk).EchoAttenuationFull));
    end % for kk

    
%
end % function          


