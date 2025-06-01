function [EchosStruct] = RecievedEchosFromSonar(FindsStruct, ObsOrPreyFlag , CurrentFreqs, nTime, CurrentPulseNum, varargin)
% [CurrentEchos] = RecivedEchosFromSonar(ObsStruct, PreyStruct):
%   This function will calculate the recived echos (Delay Times and Attenuation)  of the Sonar Pulse
% Inputs- 
%       Obsdtruct- the relative distance and angle of the findings from
%               the bat [Samples],
%       ObsOrPreyFlag = 'Prey' or 'Obs' or 'Conspecifics'- the relative distance and angle of the preys from
%               the bat
%       CurrentFreq - the current Freq transmittes
%       PulseWidth - the time length of the transmiitted Pulse [Samples]
%       xySoundVelocity - [samples]
% OutPuts-
%       CurrentEchosStrct- a Struct with the fields:
%               
%               for echo Echo - EchoSource, EchoDelayTime, EchoAttenuation,
%               EchoPulseTotalTime, 
%               TargetIndex - the inexing of the echos from input
%               
% Simplified assumptions:
%           #1  the movement of the bat is neglectable during the
%               transmitting pulse, the times are calculated from the start of the
%               pulse
%           #2  Attenuation in free Space as the RADAR Equation, Contant
%               during each pulse
%           #3  Each "point" detected in FindsStruct will Echo the full pulse 
%           #5  Prey Radar Section is const and a parmeter
%           #6  Recieving an Transmitting Gain are uniform function in the
%               BeamWidth
%           #7  Use of Multiactter Target    
% Swarm Changes by Omer Jan2021
% Updated  06Sep2022


NumberOfInputs= nargin; 
if NumberOfInputs >= 5
    AllParams =  varargin{1};
    xyResolution = AllParams.SimParams.xyResolution;
    mSoundVelocity = AllParams.SimParams.SoundV0; % Sound Velocity in meters
    TargetWingLength = AllParams.BatSonarParams.TargetWingLength;
    SampleTime = AllParams.SimParams.SampleTime;
else % if NumberOfInputs >= 5;
    xyResolution = 0.01;
    mSoundVelocity = 343;
    TargetWingLength = 0.02;
    SampleTime = 0.001;
end %if NumberOfInputs >= 5;
if NumberOfInputs == 7
    NumOfEchos =  varargin{2};
    DetectedPreys = varargin{3};
end %if NumberOfInputs == 7

% xyResolution = 0.01;
% xySoundVelocity = mSoundVelocity * xyResolution;

switch ObsOrPreyFlag
    case 'Prey'
        TargetNum  = FindsStruct.TargetsID;
        TargetArea = TargetWingLength.^2*pi; % assuming winglength ~ 1.5cm
        Distances  = FindsStruct.Distances;
        Angles     = FindsStruct.TargetAngle;
        NumOfEchos = FindsStruct.NumOfTargets;
    
    case 'Obs'
        TargetNum  = 1:numel(FindsStruct.Distances);
        TargetArea = AllParams.BatSonarParams.ClutterPointReflecetingArea; % 0.0025; %0.5e-4; % Aproximation for the area the pulse covers in samples of 0.1 radian (5cmX5cm)
        TargetArea = TargetArea * 10.^(AllParams.BatSonarParams.ClutterGainDB/10);
        Distances  = FindsStruct.Distances;
        Angles     = FindsStruct.TargetAngle;
        NumOfEchos = length(Distances);
    
    case 'Conspecifics'
        TargetWingLength = AllParams.BatSonarParams.BatWingLength; % 0.12; % winglength of the bat
        TargetArea       = TargetWingLength.^2*pi;
        % new 06Sep22
%         if strcmp(AllParams.SimParams.TestMode, 'swarm') 
%             minDist = AllParams.BatSonarParams.BatDetectionRange;
% %             minDist = AllParams.SimParams.swarm_detect_cons_dist;
%         else
%             minDist = AllParams.BatSonarParams.BatDetectionRange;
%         end
%         relevant_idx = find(FindsStruct.Distances * xyResolution <= minDist);
        relevant_idx = true(size(FindsStruct.TargetsID)); % test all distances from all conpss
        TargetNum    = FindsStruct.TargetsID(relevant_idx);
        Distances    = FindsStruct.Distances(relevant_idx);
        Angles       = FindsStruct.TargetAngle(relevant_idx);
        NumOfEchos   = numel(relevant_idx);

end %switch ObsOrPreyFlag

%%% calculations
%     NumOfEchos = length(FindsStruct.Distance); % How many differnt angles are ecohing...
    EchosStruct.NumOfEchos = NumOfEchos;
    EchosStruct.PulseDuration = length(CurrentFreqs);
    [SortedEcosDistances, IndexSorting] = sort(Distances);
    
%     EchosStruct.EchosTimes = SortedEcosDistances*2./xySoundVelocity ;   % the delays of the echos' sorted from nearest to farset
    mDistances = SortedEcosDistances*xyResolution; % Distnatnces in meters
    Angles = Angles(IndexSorting);
    EchosStruct.EchosTimes = mDistances*2./mSoundVelocity ./SampleTime ;
    MedFreq = median(CurrentFreqs); % the median frequncy for flaculation
    EchosStruct.EchosAttenuation = ...
        CalculateSignalAtten(mDistances, Angles, TargetArea, MedFreq, 'meters', mSoundVelocity, AllParams, ObsOrPreyFlag);
    
    
    EchosStruct.EchoDetailed(max(NumOfEchos,1)) = struct('TargetIndex',[], 'EchoAttenuationFull',[], 'EchosnTimesFull',[], 'EchosFreqsFull',[]);
    for kk = 1:NumOfEchos
        RoundTime = round([EchosStruct.EchosTimes(kk)]);
        EchosStruct.EchoDetailed(kk).TargetIndex = TargetNum(IndexSorting(kk));
        EchosStruct.EchoDetailed(kk).EchoAttenuationFull = ...
            CalculateSignalAtten(mDistances(kk), Angles(kk), TargetArea, CurrentFreqs, 'meters', mSoundVelocity, AllParams, ObsOrPreyFlag); 
        EchosStruct.EchoDetailed(kk).EchosnTimesFull = nTime + (0: EchosStruct.PulseDuration-1) + RoundTime;...
%             nTime + [RoundTime : ( RoundTime +EchosStruct.PulseDuration-1 )];
         EchosStruct.EchoDetailed(kk).EchosFreqsFull = CurrentFreqs;
    end % for kk

    
    EchosStruct.TransmittedPulseTime = nTime;
%     EchosStruct.EchosAngles = Angles(IndexSorting);
%     EchosStruct.EchosmDistances = mDistances;
    % update the number of the conspecifics targets to all conspecifics
    % (not only the detected ones)
    
    if strcmp(ObsOrPreyFlag, 'Conspecifics')
        EchosStruct.TargetIndex = TargetNum(IndexSorting);
    else
        EchosStruct.TargetIndex = IndexSorting; 
    end % if strcmp
    
    EchosStruct.EchosAngles = Angles;
    EchosStruct.EchosmDistances = mDistances;
    EchosStruct.EchosCenterFreq = MedFreq;
    EchosStruct.TransmittedPulseNum = CurrentPulseNum;
% % %     %%% Output
% % %     EchosStruct.NumOfEchos = NumOfEchos;
% % %     EchosStruct.EchosTimes = SortedEchosTimes;
% % %     EchosStruct.EchosAttenuation = SortedEchosAttenuation;
% % %     EchosStruct.TransmittedPulseTime = nTime;
% % %     
% % % 
end % function          
