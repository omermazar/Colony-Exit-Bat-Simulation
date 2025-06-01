function [FullFlightStages, PulsesStageStruct, PulsesStageCell, PulsesStatisticsStruct] =...
    AnalyzeFlightStages( NumOfTransmittedPulses, PulseTimesVec, ManueverPulseTypeCell , ...
    AllInterferencePulseNum, DetectionPulses, PreysDetectionMatrix , PreysInterferenceMatrix, PulseHuntedPreyVec, NumOfSamples )

%%%% Decision of the Flight Stages according to Bat Behavior:
% Inputs - 
%   NumOfTransmittedPulses
%   ManueverPulseTypeCell - cell in size (1,NumOfPulsese) Of the ManuveType for each Pulse
%   ('Foraging' / ObsMan' / 'Hunting') 
%   AllInterferencePulseNum - the PulseNumbers which were interferred
%   DetectionPulses = the pulses that preys were detectid in
%   PreysDetectionMatrix - matrix of NumOfPreyrows X NumOfPulses collums with zeros and 1ones for detections
%   PreysInterferenceMatrix - matrix od ones for interfernce
%   NumOfSamples = Total Num Of Samples in the flight
%   PulseHuntedPreyVec - the vector of the prey the bat is trying to hunt
%       in each pulse (zero if none)
% the decision rules: 
% Start with Search
% Go To search if 4 out of last 6 pulses are 'Hunting' / 'ObsMan' 
% Go Back to Search or from Hunting/ObsMan (and vv) if 2/3 is 'Foraging'


NumOfStages = 1; % total Number of stages  'Search' / 'Approach' / 'ObsManuever' 
NumOfSearchStages = 1;
SearchStagesVec = zeros(1,NumOfTransmittedPulses);
NumOfApproachStages = 0;
ApproachStagesVec = zeros(1,NumOfTransmittedPulses);
NumOfObsManueverStages = 0;
ObsManueverVec = zeros(1,NumOfTransmittedPulses);

PulsesStageCell = repmat({''},1,NumOfTransmittedPulses); % 'Search' / 'Approach' / 'ObsManuever' 
PulsesStageStruct(NumOfTransmittedPulses) = struct('StageNum',[], 'StageType',[], 'StageFirstPulse',[],'StageLastPulse',[]);

% MinRatioForSearchJamm = 0.6;

BackCounterFromSearchToApp = 7; 
BackCounterToObsMan = 1;
BackCounterFromApproachToSearch = 3;
BackCounterFromObsMan = 1;

BackCounter = 1; % In norder to avoid negative indices
% MinFindsToSwitchFromSearch = 4; 
MinFindsToSwitchFromSearchToApp = 5;
MinFindsToSwitchFromSearchToObsMan = 1; 
MinFoaragingToSwitchFromApproachToSearch = 2;
MinHuntingToSwitchFromObsMan = 1;
MinFindsToSwitchToObsMan = 1;

% NumHofHuntingInBackCounter = 0;
% NumHofObsManInBackCounter = 0;
% NumHofForagingInBackCounter = 0;

SwitchFlag =0;

%%% Init
CurrentStage = 'Search';
PulsesStageStruct(1).StageNum = 1;
PulsesStageStruct(1).StageType = CurrentStage;
PulsesStageStruct(1).StageFirstPulse = 1; 
PulsesStageStruct(1).StageLastPulse = 1;
SearchStagesVec(1) = 1;


for kPulse = 1: NumOfTransmittedPulses
 
    NumHofObsManInBackCounter = sum( strcmp('ObsMan', ManueverPulseTypeCell(kPulse-BackCounterToObsMan+1: kPulse))) ;
    
    switch CurrentStage
        case 'Search'
                % switch to ObsManuever Stage
            if NumHofObsManInBackCounter >= MinFindsToSwitchFromSearchToObsMan %
                SwitchFlag = 1;
                CurrentStage = 'ObsManuever';
                NumOfStages = NumOfStages+1;
                NumOfObsManueverStages = NumOfObsManueverStages +1;
                ObsManueverVec(NumOfObsManueverStages) = NumOfStages;
                
            else
                BackCounter = min(kPulse,BackCounterFromSearchToApp);
                LastPulsesMans = ManueverPulseTypeCell((kPulse-BackCounter+1):kPulse);
                NumHofHuntingInBackCounter = sum( strcmp('Hunting', LastPulsesMans)) ;
                
                % switch to Approach Stage
                if NumHofHuntingInBackCounter >= MinFindsToSwitchFromSearchToApp
                    SwitchFlag = 1;
                    CurrentStage = 'Approach';
                    NumOfStages = NumOfStages+1;
                    NumOfApproachStages = NumOfApproachStages +1;
                    ApproachStagesVec(NumOfApproachStages) = NumOfStages; % updating the vector of stages
                   
                    % continue Searching
                else % if NumHofHuntingInBackCounter
                    SwitchFlag = 0;
                    
                end % if NumHofHuntingInBackCounter >= MinFindsToSwitchFromSearch2App
              
            end % if NumHofHuntingInBackCounter
            
            
        case 'Approach'
            % switch to ObsManuever Stage
            if NumHofObsManInBackCounter >= MinFindsToSwitchFromSearchToObsMan %
                SwitchFlag = 1;
                CurrentStage = 'ObsManuever';
                NumOfStages = NumOfStages+1;
                NumOfObsManueverStages = NumOfObsManueverStages +1;
                ObsManueverVec(NumOfObsManueverStages) = NumOfStages;

            else
                BackCounter = min(kPulse,BackCounterFromApproachToSearch);
                LastPulsesMans = ManueverPulseTypeCell((kPulse-BackCounter+1):kPulse);
                NumHofForagingInBackCounter = sum( strcmp('Foraging', LastPulsesMans)) ;
                
                % switch to Search Stage if there are 2of3 Foraging Manuever
                if NumHofForagingInBackCounter >= MinFoaragingToSwitchFromApproachToSearch
                    SwitchFlag = 1;
                    CurrentStage = 'Search';
                    NumOfStages = NumOfStages+1;
                    NumOfSearchStages = NumOfSearchStages +1;
                    SearchStagesVec(NumOfSearchStages) = NumOfStages;
                    
                    % Keep Approaching
                else
                    SwitchFlag = 0;
                    
                end %if NumHofForagingInBackCounter >= MinFindsToSwitchFromMan
            end % if NumHofHuntingInBackCounter
            
        case 'ObsManuever'
            % Keep ObsManuevirng
            if NumHofObsManInBackCounter >= MinFindsToSwitchToObsMan %
                SwitchFlag = 0;
            
            else % if NumHofObsManInBackCounter >= MinFindsToSwitchFromSearchToObsMan %
                BackCounter = min(kPulse,BackCounterFromObsMan);
                LastPulsesMans = ManueverPulseTypeCell((kPulse-BackCounter+1):kPulse);
                NumHofHuntingInBackCounter = sum( strcmp('Hunting', LastPulsesMans)) ;
              
                % Switch to Approach
                if NumHofHuntingInBackCounter >= MinHuntingToSwitchFromObsMan
                    SwitchFlag = 1;
                    CurrentStage = 'Approach';
                    NumOfStages = NumOfStages+1;
                    NumOfApproachStages = NumOfApproachStages +1;
                    ApproachStagesVec(NumOfApproachStages) = NumOfStages; % updating the vector of stages
                    
                    % Switch to Search
                else % if NumHofHuntingInBackCounter >= MinFoaragingToSwitchFromApproachToSearch
                    SwitchFlag = 1;
                    CurrentStage = 'Search';
                    NumOfStages = NumOfStages+1;
                    NumOfSearchStages = NumOfSearchStages +1;
                    SearchStagesVec(NumOfSearchStages) = NumOfStages;
                end % if NumHofHuntingInBackCounter >= MinFoaragingToSwitchFromApproachToSearch
            end % if NumHofObsManInBackCounter >= MinFindsToSwitchFromSearchToObsMan %
                            
    end %switch CurrentStage
    
    % outputs
    PulsesStageCell{kPulse} = CurrentStage;
    
    if SwitchFlag
        PulsesStageStruct(NumOfStages).StageNum = NumOfStages;
        PulsesStageStruct(NumOfStages).StageType = CurrentStage;
        PulsesStageStruct(NumOfStages).StageFirstPulse = kPulse;
        PulsesStageStruct(NumOfStages).StageLastPulse = kPulse;
        
    else % if SwitchFlag
        PulsesStageStruct(NumOfStages).StageLastPulse = kPulse;
        
    end % if SwitchFlag
    
    if BackCounter < BackCounterFromSearchToApp
        BackCounter = BackCounter+1;
    end % if BackCounter == NumHofHuntingInBackCounter
    
end % for kk = 1: NumOfTransmittedPulses

%%% Anylizing interference ratio in each Stage
NumberOfSuccesulSearchStages = 0;
NumberOfJammedSearchStages = 0;
NumberOfSuccesulApproachStages = 0;
NumberOfJammedApproachStages = 0;

MaxSequenceJamsForSearch = 3;
MaxSequenceJamsForApproach = 3;

MinDetectionForSuccesfulSearch = 4;
MinRatioForApproachJamm = 0.33;
MinDetectionForSuccesfulApproach = 5;

for nn= 1:NumOfStages
    minPulse = PulsesStageStruct(nn).StageFirstPulse;
    maxPulse = PulsesStageStruct(nn).StageLastPulse;
    PulsesStageStruct(nn).StageStartTIme = PulseTimesVec(minPulse);
    PulsesStageStruct(nn).StageEndTime = PulseTimesVec(maxPulse)-1;
    InterferencePulseIndex = find((AllInterferencePulseNum >= minPulse) & (AllInterferencePulseNum <= maxPulse));
    PulsesStageStruct(nn).NumOfInterInStage = length(InterferencePulseIndex);
    PulsesStageStruct(nn).NumOfDetections = length(find((DetectionPulses >= minPulse) & (DetectionPulses <= maxPulse)));
    
    if PulsesStageStruct(nn).NumOfDetections >0
        PulsesStageStruct(nn).InterferenseRatio = PulsesStageStruct(nn).NumOfInterInStage ./ PulsesStageStruct(nn).NumOfDetections;
    else % if 
         PulsesStageStruct(nn).InterferenseRatio = [];
    end % if NumOfdetections
    
   
    %%% Analyze the suuccess and jamming of each stage %%%
    
    PulsesStageStruct(nn).JammedStageFlag = 0;
    PulsesStageStruct(nn).SuccessfulStageFlag = 0;
    PreysDetectedMatrixStage = PreysDetectionMatrix(:, minPulse:maxPulse );
    PreysInterferenceMatrixStage = PreysInterferenceMatrix(:, minPulse:maxPulse ); 
    HuntedPreyInStage = PulseHuntedPreyVec( minPulse : maxPulse );

    
%     Is2SequenceOfJamFlag = ~isempty(find(diff(InterferencePulseIndex) == 1,1));
    if PulsesStageStruct(nn).NumOfDetections >0
        switch PulsesStageStruct(nn).StageType
            case 'Approach'
                %%% Success Rule- at least 5 detections of hunted prey (before jamming)
                %%% Jamming- if more tha 33% of detections of the hunted
                %%% are interferres  or 2 consecutive detection are lost
                
                %%% Success Condition
                UniqueHuntedPreyDetected = unique(nonzeros(HuntedPreyInStage));
                NumOfDetectionOfEachPrey = nonzeros(accumarray( nonzeros(HuntedPreyInStage') ,1))';
                IndexOfSuccess = find (NumOfDetectionOfEachPrey >= MinDetectionForSuccesfulApproach);
                
                SuccessfulHuntedPreysDetected = UniqueHuntedPreyDetected(IndexOfSuccess);
                NumOfDetectionsOfSuccessfulHuntedPreys = NumOfDetectionOfEachPrey(IndexOfSuccess);
                NumOfSucessfulHuntedPreys = length(IndexOfSuccess);
                
                if NumOfSucessfulHuntedPreys >= 1
                    PulsesStageStruct(nn).SuccessfulStageFlag = 1;
                    NumberOfSuccesulApproachStages = NumberOfSuccesulApproachStages + 1; %
                end % if max(  NumOfDetectionOfEachPrey)
                
                %%% Jamming Condition
                if PulsesStageStruct(nn).SuccessfulStageFlag >= 1
                    
                    HuntedPreysDetectionMatrix = zeros(length(SuccessfulHuntedPreysDetected), length(HuntedPreyInStage));
                    HuntedPreysInterfernceMatrix = HuntedPreysDetectionMatrix;
                    for jj = 1:length(SuccessfulHuntedPreysDetected)
                        IndexA = find(HuntedPreyInStage == SuccessfulHuntedPreysDetected(jj));
                        HuntedPreysDetectionMatrix(jj, IndexA) = ...
                            PreysDetectedMatrixStage(SuccessfulHuntedPreysDetected(jj) ,IndexA );
                        HuntedPreysInterfernceMatrix(jj, IndexA) = ...
                            PreysInterferenceMatrixStage(SuccessfulHuntedPreysDetected(jj) , IndexA );
                    end % for jj
                    NotInterferredPreysMatrix = HuntedPreysDetectionMatrix - HuntedPreysInterfernceMatrix;
                    NumOfNotInterferredHuntedPrey = sum(NotInterferredPreysMatrix');
                    ApproacInterfernceRatio = (NumOfDetectionsOfSuccessfulHuntedPreys - NumOfNotInterferredHuntedPrey) ./ NumOfDetectionsOfSuccessfulHuntedPreys;
                    JammedHuntedPreys = find(ApproacInterfernceRatio >= MinRatioForApproachJamm);
                    % Cond1
                    if length(JammedHuntedPreys) == NumOfSucessfulHuntedPreys
                       PulsesStageStruct(nn).JammedStageFlag = 1;
                       NumberOfJammedApproachStages = NumberOfJammedApproachStages+1;
                    
                    % Cond2
                    else % if max (ApproacInterfernceRatio)
                        
                        IndexUnInterferredPreys = find(NumOfNotInterferredHuntedPrey >= MinDetectionForSuccesfulApproach);
                        UnInterferredPreys = SuccessfulHuntedPreysDetected(IndexUnInterferredPreys);
                        JamMat = zeros(length(UnInterferredPreys), length(HuntedPreyInStage));
                        for jj = 1:length(UnInterferredPreys)
                            IndexB = find(HuntedPreyInStage == UnInterferredPreys(jj));
                            JamMat(jj,IndexB) = HuntedPreysInterfernceMatrix(IndexUnInterferredPreys(jj), IndexB);
                        end % for jj
                        NumOfUnInterferred = CheckSequnceOfInterference(UnInterferredPreys, JamMat, MaxSequenceJamsForApproach);
                        
                        if NumOfUnInterferred < 1
                           PulsesStageStruct(nn).JammedStageFlag = 1;
                           NumberOfJammedSearchStages = NumberOfJammedSearchStages+1;
                        end %if NumOfUnInterferred < 1
                        
                    end % if max (ApproacInterfernceRatio)
                    
                end % if ulsesStageStruct(nn).SuccessfulStageFlag
%                 
            case 'Search'
                %%% Success Rule (before jamming)  - at least 4 detection
                %%% of the same prey in the STAGE
                %%% Jamming Rules (after jamming) - 
                %%% 1) less than 4 detections
                %%% are not jammed - or- 
                %%% 2) at least 2 consequative echos of
                %%% the same prey
                
                 %%% Success condition
                SumOfPreysDetectionsInStageVec = sum(PreysDetectedMatrixStage');
                MaxDetectionOfPrey = max(SumOfPreysDetectionsInStageVec);
                if MaxDetectionOfPrey >= MinDetectionForSuccesfulSearch
                    PulsesStageStruct(nn).SuccessfulStageFlag = 1;
                    NumberOfSuccesulSearchStages = NumberOfSuccesulSearchStages+1;
                end% if MaxDetectionOfPrey >= MinFindsToSwitchFromSearchToApp
                
                %%% Jamming conditions
                if PulsesStageStruct(nn).SuccessfulStageFlag 
                    NotInterferredMatrix = PreysDetectedMatrixStage - PreysInterferenceMatrixStage;
                    SumOfNotInterferred = sum(NotInterferredMatrix');
                % Cond1 
                    if max(SumOfNotInterferred) < MinDetectionForSuccesfulSearch
                       PulsesStageStruct(nn).JammedStageFlag = 1; 
                       NumberOfJammedSearchStages = NumberOfJammedSearchStages+1;
                % Cond2
                    else % max(SumOfNotInterferred)
                        
                        UnInterferredPreys = find(SumOfNotInterferred >= MinDetectionForSuccesfulSearch);
                        JamMat = PreysInterferenceMatrixStage(UnInterferredPreys,:);
                        
                        NumOfUnInterferred = CheckSequnceOfInterference(UnInterferredPreys, JamMat, MaxSequenceJamsForSearch);
                        
                        if NumOfUnInterferred < 1
                           PulsesStageStruct(nn).JammedStageFlag = 1;
                           NumberOfJammedSearchStages = NumberOfJammedSearchStages+1;
                        end %if NumOfUnInterferred < 1
                        
                    end % if max(SumOfPreysDetectionsInStageVec) < MinFindsToSwitchFromSearchToApp
                    
                    
                end%if PulsesStageStruct(nn).SuccessfulStageFlag 
                
                
        end % switch PulsesStageStruct(NumOfStages).StageType(nn)
        
    end % if PulsesStageStruct(nn).NumOfDetections >0
    
end% for nn= 1:NumOfStages


%%% OutPuts

FullFlightStages = Pulses2Flight(PulsesStageCell, PulseTimesVec, NumOfSamples);

PulsesStatisticsStruct.NumOfStages= NumOfStages;
PulsesStatisticsStruct.NumOfSearchStages= NumOfSearchStages;
PulsesStatisticsStruct.NumOfObsManueverStages= NumOfObsManueverStages;
PulsesStatisticsStruct.NumOfApproachStages= NumOfApproachStages;

PulsesStatisticsStruct.SearchStagesVec = int16(nonzeros(SearchStagesVec)) ;
PulsesStatisticsStruct.ApproachStagesVec = int16(nonzeros(ApproachStagesVec)) ;
PulsesStatisticsStruct.ObsManueverVec = int16(nonzeros(ObsManueverVec)) ;

PulsesStatisticsStruct.NumberOfSuccesulSearchStages = NumberOfSuccesulSearchStages;
PulsesStatisticsStruct.NumberOfJammedSearchStages = NumberOfJammedSearchStages;
if NumberOfSuccesulSearchStages
    PulsesStatisticsStruct.JammedToSuccessfulSearchStages = NumberOfJammedSearchStages/ NumberOfSuccesulSearchStages;
else
    PulsesStatisticsStruct.JammedToSuccessfulSearchStages =[];
end % if NumberOfSuccesulSearchStages

PulsesStatisticsStruct.NumberOfSuccesulApproachStages = NumberOfSuccesulApproachStages;
PulsesStatisticsStruct.NumberOfJammedApproachStages = NumberOfJammedApproachStages;
if NumberOfSuccesulApproachStages
    PulsesStatisticsStruct.JammedToSuccessfulApproachStages = NumberOfJammedApproachStages/ NumberOfSuccesulApproachStages;
else
    PulsesStatisticsStruct.JammedToSuccessfulApproachStages =[];
end % if NumberOfSuccesulSearchStages

    
end % main function



%%%%%%%%%%%%% 

% 
% function [FullVec] = Pulses2Flight(PulsesStageCell, StartPulseTimes, NumOfSamples)
% %%% returns the stage in everysample of the bat from the stage in each pulse
% 
%     FullVec= repmat({''},1,NumOfSamples);
% 
%     % FullVec(StartPulseTimes) = FlightStageCell;
%     for kk = 1:length(StartPulseTimes)-1
%         FullVec(StartPulseTimes(kk): StartPulseTimes(kk+1)-1 ) = PulsesStageCell(kk);
%     end % for kk
%     FullVec(StartPulseTimes(end):end) = PulsesStageCell(end);
% end % function Pulses2Flight


%%%%%%%%%%%%%%%%%%


function [NumOfUnInterferred] = CheckSequnceOfInterference(UnInterferredPreys, JamMat, MaxSequlnceLength)
%%% Returns the length of sequnces of interfernces in the stage

    NumOfUnInterferred = length(UnInterferredPreys);

    for kk = 1:length(UnInterferredPreys)
        JamIndex = find(JamMat(kk,:));
        DiffIndex = diff(JamIndex);
        B = find([DiffIndex, inf] > 1);
        JammingSequencesLength = diff([0,B]);
        if max(JammingSequencesLength) >= MaxSequlnceLength % there are at least two consequent jamming
            NumOfUnInterferred = NumOfUnInterferred -1;
        end % if ~isempty(ConsequncesInterfernce)
    end % for kk

end % function