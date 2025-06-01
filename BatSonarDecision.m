function [BAT] = BatSonarDecision(BAT, CurrPulseNum, AllParams)
% This fuction will decide the PRF (IPI), the Pulse Duration and he Pulse Power of the
% Sonar, According to the manuever command decicded
% Inputs: BAT.ManueverCmdStruct
%         Current PulseNum
%          All Sonar Parmaters
%         BAT.CurrentTimeToSonar(nTime) and BAT.CurrentPulseWidth(nTime)
%
% Outputs:  CurrentTimeToSonar(nTime+1) - (Current PRF in n Sample Time)
%           CurrentPRF
%           CurrentPulseWidth(nTime+1) - not integer
%           CurrentPulseDuration = round(CurrentPulseWidth)
%           CurrentChirp
%           ChirpStep
%           PulsePower

%%% Changed on APR-2022 for Moth Jam
% Control on the Terminal Freq too during the Appraoch and Buzz


% General Parameters
SampleTime      = AllParams.SimParams.SampleTime;
xyResolution    = AllParams.SimParams.xyResolution;
xySoundVelocity = AllParams.SimParams.SoundV0 * SampleTime/xyResolution;
NumberOfSteps   = 10;

%%% BAT SONAR Paramaters
%%
DetectionRange     = AllParams.BatSonarParams.BatDetectionRange/xyResolution;
SensorAnalyzeDelay = 1e-3/SampleTime;

% PulsePower
PulseNominalPower = AllParams.BatSonarParams.PulsePower;
PulseMinPower     = AllParams.BatSonarParams.PulseMinPower;

% Terminal Frequency - NEW Apr2022 - Moth JAM
TerminalFreq_Search        = AllParams.BatSonarParams.TerminalFreq_Search        + BAT.TerminalFreqUnique;
TerminalFreq_ApproachStart = AllParams.BatSonarParams.TerminalFreq_ApproachStart + BAT.TerminalFreqUnique;
TerminalFreq_ApproachEnd   = AllParams.BatSonarParams.TerminalFreq_ApproachEnd   + BAT.TerminalFreqUnique;
TerminalFreq_BuzzStart     = AllParams.BatSonarParams.TerminalFreq_BuzzStart     + BAT.TerminalFreqUnique;
TerminalFreq_BuzzEnd       = AllParams.BatSonarParams.TerminalFreq_BuzzEnd       + BAT.TerminalFreqUnique;
TerminalFreq_BuzzFinal     = AllParams.BatSonarParams.TerminalFreq_BuzzFinal     + BAT.TerminalFreqFinalBuzzUnique;

% IPI Paramaeters
IPI_Search        = AllParams.BatSonarParams.IPI_Search *1e-3/SampleTime;
IPI_ApproachStart = AllParams.BatSonarParams.IPI_ApproachStart *1e-3/SampleTime;
IPI_ApproachEnd   = AllParams.BatSonarParams.IPI_ApproachEnd *1e-3/SampleTime;
IPI_BuzzStart     = AllParams.BatSonarParams.IPI_BuzzStart *1e-3/SampleTime;
IPI_BuzzEnd       = AllParams.BatSonarParams.IPI_BuzzEnd *1e-3/SampleTime;
IPI_BuzzFinal     = AllParams.BatSonarParams.IPI_BuzzFinal *1e-3/SampleTime;

% Pulse Duration
PulseDuration_Search       = AllParams.BatSonarParams.PulseDuration_Search *1e-3/SampleTime;
PulseDuration_ApprachStart = AllParams.BatSonarParams.PulseDuration_ApprachStart *1e-3/SampleTime;
PulseDuration_ApproachEnd  = AllParams.BatSonarParams.PulseDuration_ApproachEnd *1e-3/SampleTime;
PulseDuration_BuzzStart    = AllParams.BatSonarParams.PulseDuration_BuzzStart *1e-3/SampleTime;
PulseDuration_BuzzEnd      = AllParams.BatSonarParams.PulseDuration_BuzzEnd *1e-3/SampleTime;
PulseDuration_BuzzFinal    = AllParams.BatSonarParams.PulseDuration_BuzzFinal *1e-3/SampleTime;

% PulseBandwidth kHz
ChirpSpan_Search        = AllParams.BatSonarParams.ChirpSpan_Search ;
ChirpSpan_ApproachStart = AllParams.BatSonarParams.ChirpSpan_ApproachStart;
ChirpSpan_ApproachEnd   = AllParams.BatSonarParams.ChirpSpan_ApproachEnd;
ChirpSpan_Buzz          = AllParams.BatSonarParams.ChirpSpan_Buzz;
ChirpSpan_BuzzFinal     = AllParams.BatSonarParams.ChirpSpan_BuzzFinal;

% obstacles Behavior dervied from hunt
MinObsTimeToSonar = IPI_ApproachEnd;
ShortObsPulse     = PulseDuration_BuzzEnd;
StepPRF           = (IPI_Search- IPI_ApproachEnd)./NumberOfSteps;
StepPulseTime     =  (PulseDuration_Search- PulseDuration_ApproachEnd)./NumberOfSteps;
StepPulsePower    = (PulseNominalPower - PulseMinPower)./NumberOfSteps;

PulseRecievedTimeToAnalyze = 2*DetectionRange/xySoundVelocity + SensorAnalyzeDelay;

% random Search_IPI
if sum(strcmp('RandomIPIFlag', fields(AllParams.BatSonarParams)))
    if rand(1,1) < AllParams.BatSonarParams.RandomIPIChanceofShort
        PulseDuration_Search = round(PulseDuration_Search/2);
        IPI_Search = round(IPI_Search/2);
    end % if rand
end % if sum(strcmp())

%%% OldModel
% MinTimeToSonar = 1/AllParams.BatSonarParams.BatMaxPRF/SampleTime;
% RegTimeToSonar = 1./ AllParams.BatSonarParams.BatNominalPRF /SampleTime; % regular- in free terrain
% MinObsTimeToSonar = MinTimeToSonar *2; %%
% LongPulse = AllParams.BatSonarParams.PulseTimeLong / SampleTime;
% ShortPulse = AllParams.BatSonarParams.PulseTimeShort/SampleTime;
% ShortObsPulse = ShortPulse*2;
% StepPRF = (RegTimeToSonar- MinTimeToSonar)./NumberOfSteps;
% StepPulseTime=  (LongPulse- ShortPulse)./NumberOfSteps;
% PulseNominalPower = AllParams.BatSonarParams.PulsePower;
% PulseMinPower = AllParams.BatSonarParams.PulseMinPower;
%
%
% StepPulsePower = (PulseNominalPower - PulseMinPower)./NumberOfSteps;
%
% PulseRecievedTimeToAnalyze = 2*DetectionRange/xySoundVelocity + SensorAnalyzeDelay;
% ChirpFlag = AllParams.BatSonarParams.ChirpFlag;
% ChirpSpan = AllParams.BatSonarParams.ChirpSpan*ChirpFlag; % zero if not chirping
% ChirpLowestFreq = AllParams.BatSonarParams.CenterFreq - ChirpSpan/2;

%%

%%% update AllPARAMS
% try
DistanceToApproach = AllParams.BatFlightParams.DistanceFromPreyToApproach./xyResolution;
DistanceToBuzz =  AllParams.BatFlightParams.DistanceToBuzz./xyResolution;
DistanceToFinalBuzz = AllParams.BatFlightParams.DistanceToFinalBuzz/xyResolution;

%%%% NEW June22 - caveExit update
if strcmp(AllParams.SimParams.TestMode, 'caveExit')
    DistanceToApproach = AllParams.BatSonarParams.BatDetectionRange ./ xyResolution + 1;
end
% catch
%     DistanceFromPreyToChangeSonar = 2./xyResolution;
%     DistanceToBuzz = 0.4./xyResolution;
% end % try
DiffX = DistanceToApproach - DistanceToBuzz;

%%% Input Parameters
% BAT.ManueverType;  %'Hunting' or 'ObsMan' or 'Foraging'
%BAT.ManueverPower;

% ACTION
BAT.PulseFreqsCommand = -1;
FinalBuzzFlag = 0;

% TerminalFreq = AllParams.BatSonarParams.TerminalFreq ; % the lowest feq in the chirp [kHz]
% TerminalFreq = BAT.TerminalFreqUnique;
TerminalFreq = TerminalFreq_Search;

%%% XXXXXXX %%%%%%
ManueverCmdStruct = BAT.ManueverCmdStruct(CurrPulseNum);

switch ManueverCmdStruct.ManueverType
    case {'Hunting', 'ExitMan', 'FollowWall', 'FollowWallCarefully'}
        Dist2HuntedPrey = ManueverCmdStruct.Dist2HuntedPrey;

        switch ManueverCmdStruct.ManueverStage % 'Search' / 'Approach' /'Buzz' / 'ObstacleManuver'/ 'AvoidBatMan'
            % sonar paramters as function of distance from prey
            case 'Search'
                BAT.CurrentTimeToSonar = IPI_Search;
                BAT.CurrentPulseWidth = PulseDuration_Search;
                BAT.PulsePower = PulseNominalPower;
                FreqBandWidth = ChirpSpan_Search;

            case {'Approach', 'CaveExit'}
                Dist2HuntedPrey = min(Dist2HuntedPrey,DistanceToApproach);
                TerminalFreq = interp1( [DistanceToApproach, DistanceToBuzz], ...
                    [TerminalFreq_ApproachStart, TerminalFreq_ApproachEnd], Dist2HuntedPrey, 'linear');

                BAT.CurrentTimeToSonar = interp1( [DistanceToApproach, DistanceToBuzz], ...
                    [IPI_ApproachStart, IPI_ApproachEnd], Dist2HuntedPrey, 'linear');

                BAT.CurrentPulseWidth = interp1( [DistanceToApproach, DistanceToBuzz], ...
                    [PulseDuration_ApprachStart, PulseDuration_ApproachEnd], Dist2HuntedPrey, 'linear');

                BAT.PulsePower = interp1( [DistanceToApproach, DistanceToBuzz], ...
                    [PulseNominalPower, PulseMinPower], Dist2HuntedPrey, 'linear');

                FreqBandWidth = interp1( [DistanceToApproach, DistanceToBuzz], ...
                    [ChirpSpan_ApproachStart, ChirpSpan_ApproachEnd], Dist2HuntedPrey, 'linear');

                %                 BAT.CurrentTimeToSonar = MinTimeToSonar + ...
                %                     max(0, (RegTimeToSonar - MinTimeToSonar)/DiffX * (Dist2HuntedPrey - DistanceToBuzz) ) ;
                %                 BAT.CurrentPulseWidth = ShortPulse + ...
                %                     max(0, (LongPulse - ShortPulse)/DiffX * (Dist2HuntedPrey - DistanceToBuzz) );
                %                 BAT.PulsePower = PulseMinPower + ...
                %                     max(0, (PulseNominalPower - PulseMinPower)/DiffX * (Dist2HuntedPrey - DistanceToBuzz) );

            case 'Buzz'
                %                 BAT.PulsePower  =  interp1( [ DistanceToBuzz , DistanceToFinalBuzz], ...
                %                         [PulseMinPower, PulseMinPower-10], Dist2HuntedPrey, 'linear');

                if Dist2HuntedPrey >= DistanceToApproach

                    TerminalFreq = TerminalFreq_ApproachStart;
                    BAT.PulsePower  =  PulseNominalPower;
                    BAT.CurrentTimeToSonar = IPI_ApproachStart;
                    BAT.CurrentPulseWidth = PulseDuration_ApprachStart;
                    FreqBandWidth = ChirpSpan_ApproachStart;

                elseif Dist2HuntedPrey >= DistanceToBuzz
                    TerminalFreq = interp1( [DistanceToApproach, DistanceToBuzz], ...
                        [TerminalFreq_ApproachStart, TerminalFreq_ApproachEnd], Dist2HuntedPrey, 'linear');

                    BAT.CurrentTimeToSonar = interp1( [DistanceToApproach, DistanceToBuzz], ...
                        [IPI_ApproachStart, IPI_ApproachEnd], Dist2HuntedPrey, 'linear');

                    BAT.CurrentPulseWidth = interp1( [DistanceToApproach, DistanceToBuzz], ...
                        [PulseDuration_ApprachStart, PulseDuration_ApproachEnd], Dist2HuntedPrey, 'linear');

                    BAT.PulsePower = interp1( [DistanceToApproach, DistanceToBuzz], ...
                        [PulseNominalPower, PulseMinPower], Dist2HuntedPrey, 'linear');

                    FreqBandWidth = interp1( [DistanceToApproach, DistanceToBuzz], ...
                        [ChirpSpan_ApproachStart, ChirpSpan_ApproachEnd], Dist2HuntedPrey, 'linear');
                elseif Dist2HuntedPrey >= DistanceToFinalBuzz
                    TerminalFreq = interp1( [DistanceToApproach, DistanceToFinalBuzz], ...
                        [TerminalFreq_BuzzStart, TerminalFreq_BuzzEnd], Dist2HuntedPrey, 'linear');

                    BAT.PulsePower  =  interp1( [ DistanceToBuzz , DistanceToFinalBuzz], ...
                        [PulseMinPower, PulseMinPower-10], Dist2HuntedPrey, 'linear');

                    BAT.CurrentTimeToSonar = interp1( [ DistanceToBuzz , DistanceToFinalBuzz], ...
                        [IPI_BuzzStart, IPI_BuzzEnd], Dist2HuntedPrey, 'linear');

                    BAT.CurrentPulseWidth = interp1( [ DistanceToBuzz , DistanceToFinalBuzz], ...
                        [PulseDuration_BuzzStart, PulseDuration_BuzzEnd], Dist2HuntedPrey, 'linear');

                    FreqBandWidth = ChirpSpan_Buzz;

                    % FinalBuzz
                else % if Dist2HuntedPrey >= DistanceToFinalBuzz
                    BAT.PulsePower         =  PulseMinPower-10;
                    BAT.CurrentTimeToSonar = IPI_BuzzFinal;
                    BAT.CurrentPulseWidth  = PulseDuration_BuzzFinal;
                    FreqBandWidth          = ChirpSpan_BuzzFinal;
                    TerminalFreq           = TerminalFreq_BuzzFinal; % BAT.TerminalFreqFinalBuzzUnique;
                    FinalBuzzFlag          = 1;

                end %if Dist2HuntedPrey >= DistanceToFinalBuzz

                %%% DEBUGGING
                if isnan(BAT.PulsePower)
                    HOP = 'BATSONar BUZZ NAN PULSE POWER'
                end

                %                 BAT.CurrentTimeToSonar = ceil(MinTimeToSonar/4);
                %                 BAT.CurrentPulseWidth = ceil(ShortPulse/3);
                %                 BAT.PulsePower  = PulseMinPower;

        end % switch ManueverCmdStruct.ManueverStage

    case {'ObsMan', 'ObsManByMasking'} % switch ManueverCmdStruct.ManueverType
        BAT.PulsePower = max(BAT.PulsePower - StepPulsePower, PulseMinPower);

        switch  ManueverCmdStruct.ManueverPower   % 'RegularManuever' or 'CrushAvoidance'
            case 'RegularManuever'
                BAT.CurrentTimeToSonar = max(MinObsTimeToSonar, BAT.CurrentTimeToSonar - StepPRF);
                BAT.CurrentPulseWidth = max(BAT.CurrentPulseWidth-StepPulseTime, ShortObsPulse);
                FreqBandWidth = ChirpSpan_ApproachStart;
            case 'CrushAvoidance'
                % cHANgE nEW June2024
                BAT.CurrentTimeToSonar = IPI_ApproachStart; % IPI_ApproachEnd;
                BAT.CurrentPulseWidth =  PulseDuration_ApprachStart; % PulseDuration_ApproachEnd;
                BAT.PulsePower = PulseNominalPower; % PulseMinPower;
                FreqBandWidth = ChirpSpan_ApproachStart; %ChirpSpan_Buzz;
        end % switch  ManueverCmdStruct.ManueverPower

    case 'AvoidBatMan'  % switch ManueverCmdStruct.ManueverType
        if  strcmp(ManueverCmdStruct.ManueverPower, 'Buzz')
            BAT.CurrentTimeToSonar = IPI_BuzzEnd;
            BAT.CurrentPulseWidth = PulseDuration_BuzzEnd;
            BAT.PulsePower  = PulseMinPower;
            FreqBandWidth = ChirpSpan_Buzz;
        else
            Dist2Bat = min(ManueverCmdStruct.Dist2Bat, DistanceToApproach);
            TerminalFreq = interp1( [DistanceToApproach, DistanceToBuzz], ...
                [TerminalFreq_ApproachStart, TerminalFreq_ApproachEnd], Dist2Bat, 'linear');

            BAT.CurrentTimeToSonar = interp1( [DistanceToApproach, DistanceToBuzz], ...
                [IPI_ApproachStart, IPI_ApproachEnd], Dist2Bat, 'linear');

            BAT.CurrentPulseWidth = interp1( [DistanceToApproach, DistanceToBuzz], ...
                [PulseDuration_ApprachStart, PulseDuration_ApproachEnd], Dist2Bat, 'linear');

            BAT.PulsePower = interp1( [DistanceToApproach, DistanceToBuzz], ...
                [PulseNominalPower, PulseMinPower], Dist2Bat, 'linear');

            FreqBandWidth = interp1( [DistanceToApproach, DistanceToBuzz], ...
                [ChirpSpan_ApproachStart, ChirpSpan_ApproachEnd], Dist2Bat, 'linear');

        end %

        %         BAT.CurrentTimeToSonar = ceil(MinTimeToSonar/4);
        %         BAT.CurrentPulseWidth = ceil(ShortPulse/3);
        %         BAT.PulsePower  = PulseMinPower;

    case 'Foraging' % switch ManueverCmdStruct.ManueverType
        BAT.CurrentTimeToSonar = IPI_Search;
        BAT.CurrentPulseWidth = PulseDuration_Search;
        BAT.PulsePower = PulseNominalPower;
        FreqBandWidth = ChirpSpan_Search;

        %         BAT.CurrentTimeToSonar = RegTimeToSonar;
        %         BAT.CurrentPulseWidth = LongPulse;
        %         BAT.PulsePower = PulseNominalPower;

end % switch BAT.ManueverType

%% react to other bats
% chanee echo behavior as function of distance from the bat (like behavior in aparoach
% phase)

if ~strcmp(ManueverCmdStruct.ManueverStage, 'Buzz')&& ...
        ~strcmp(ManueverCmdStruct.ManueverType,'ObsMan') && ...
        ~strcmp(ManueverCmdStruct.ManueverType,'AvoidBatMan')
    if ManueverCmdStruct.ReactToBat
        Dist2Calc = max( ManueverCmdStruct.Dist2ReactBat, DistanceToBuzz);
        BAT.CurrentTimeToSonar = interp1( [DistanceToApproach, DistanceToBuzz], ...
            [IPI_ApproachStart, IPI_ApproachEnd],Dist2Calc, 'linear');

        BAT.CurrentPulseWidth = interp1( [DistanceToApproach, DistanceToBuzz], ...
            [PulseDuration_ApprachStart, PulseDuration_ApproachEnd],  Dist2Calc, 'linear');

        BAT.PulsePower = interp1( [DistanceToApproach, DistanceToBuzz], ...
            [PulseNominalPower, PulseMinPower],  Dist2Calc, 'linear');

        FreqBandWidth = interp1( [DistanceToApproach, DistanceToBuzz], ...
            [ChirpSpan_ApproachStart, ChirpSpan_ApproachEnd],  Dist2Calc, 'linear');

    end  % if ManueverCmdStruct.ReactToBat
end % if ~strcmp ...
%%

% PulseWidth must be an integer for samples
BAT.CurrentPulseWidth = max(1, BAT.CurrentPulseWidth); % cannot be less than 1
BAT.CurrentPulseDuration = round(BAT.CurrentPulseWidth);
% Decision before next sonar
BAT.CurrentTimeToSonar = round( max(IPI_BuzzFinal, BAT.CurrentTimeToSonar)); % cannot be lower than the minumu allowed

FreqBandWidth = max(1,FreqBandWidth);
%%% Updating Running Parameters %%%

%%% JAR Behavior
%check if JAR mod is on	m Start  on 2nd pulse
if (AllParams.BatSonarParams.JARMode) && (CurrPulseNum > 1)

    CurrTerminalFreq= BAT.TransmittedPulsesStruct(CurrPulseNum).ChirpMinMaxFreqs(1);
    % check if hunted prey jammed
    if BAT.PreyFindsStruct(CurrPulseNum).IsHuntedPreyMasked
        % set the JAR counter
        BAT.JARPulsesLeft = AllParams.BatSonarParams.JARMem;

        % Check Boundaries
        if FinalBuzzFlag
            MaxTerFreq = BAT.TerminalFreqFinalBuzzUnique + ...
                AllParams.BatSonarParams.TerminalFreqVariance;
            MinTerFreq = BAT.TerminalFreqFinalBuzzUnique - ...
                AllParams.BatSonarParams.TerminalFreqVariance;
        else % if FinalBuzzFlag
            MaxTerFreq = BAT.TerminalFreqUnique +...
                AllParams.BatSonarParams.TerminalFreqVariance;
            MinTerFreq = BAT.TerminalFreqUnique -...
                AllParams.BatSonarParams.TerminalFreqVariance;
        end % if FinalBuzzFlag

        % Check whether to raise or decrease the terminal freq

        % diffrence betweeen correlation and filterba
        if numel(BAT.FindsMaskingStruct(CurrPulseNum).InterMinFreq) > 1
            first_mask=  BAT.FindsMaskingStruct(CurrPulseNum).masked_prey_idx(1);
            InterMinFreq = BAT.FindsMaskingStruct(CurrPulseNum).InterMinFreq(first_mask);
        else
            InterMinFreq = BAT.FindsMaskingStruct(CurrPulseNum).InterMinFreq;
        end %numel(BAT.FindsMaskingStruct(CurrPulseNum).InterMinFreq) > 1

        if InterMinFreq <=  CurrTerminalFreq
            TerminalFreq = ...
                min(CurrTerminalFreq + AllParams.BatSonarParams.JARDelta, MaxTerFreq);
        else % if BAT.FindsMaskingStruct(CurrPulseNum).InterMinFreq <=  CurrTerminalFreq
            TerminalFreq = ...
                max(CurrTerminalFreq - AllParams.BatSonarParams.JARDelta, MinTerFreq);
        end % if BAT.FindsMaskingStruct(CurrPulseNum).InterMinFreq <=  CurrTerminalFreq



    else % if BAT.PreyFindsStruct(CurrPulseNum).IsHuntedPreyMasked
        % check if bat is allready in JAR Active
        if BAT.JARPulsesLeft > 0
            BAT.JARPulsesLeft = max(0, BAT.JARPulsesLeft - 1); % update counter
            TerminalFreq = CurrTerminalFreq;
        else % if BAT.JARPulsesLeft > 0
            % do nothing , Terminal freq is allready set
        end % if BAT.JARPulsesLeft >0
    end % %    if BAT.PreyFindsStruct(CurrPulseNum).IsHuntedPreyMasked
end % if AllParams.BatSonarParams.JARMode

[BAT.PulseFreqsCommand, BAT.ChirpMinMaxFreqs] = ...
    CalcFreqsForPulse(BAT.CurrentPulseDuration, TerminalFreq, FreqBandWidth, 'Logarithmic'); %%%% XXXXX %%%


% catch
%     pop = 'fuck!! CalcFreqs'
% end

% Decision before next sonar
BAT.PulseRecievedTimeToAnalyze =  BAT.CurrentTimeToSonar-1;

% % THis is the correct line


% BAT.ChirpStep = FreqBandWidth/(BAT.CurrentPulseWidth-1);
