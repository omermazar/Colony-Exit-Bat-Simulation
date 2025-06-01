function [SonogramMatrix, SonarEchosMat, PreyEchosVec, ObsEchosVec, InterferenceSonogram ] = ...
    BuildSonogram( NumOfTimesSonarTransmits, TransmittedPulsesStruct ...
    ,NumOfTimesObsFinds ,EchosFromObsStruct , ...
    NumOfTimesPreyFinds, EchosFromPreyStruct ,PulseWidthVec,TransPulsePower, AllParams)

    
    % Init
    SimulationTime = AllParams.SimParams.SimulationTime;
    SampleTime = AllParams.SimParams.SampleTime;
    NumOfSamples = SimulationTime/SampleTime;
    MaxFreq = AllParams.BatSonarParams.TerminalFreq + ...
        AllParams.BatSonarParams.ChirpSpan_ApproachStart + 10;
%     PulsePower = 10.^(TransPulsePower/10);
    FreqRes = 1;
    FreqGrid = unique(round([1:MaxFreq+10]/FreqRes))*FreqRes;
%     NoiseLevel = 0.1; % [SPL - linear]
    NoiseLevel = 10.^(AllParams.BatSonarParams.NoiseLeveldB/10);

    SonogramMatrix = NoiseLevel.*ones(length(FreqGrid) ,NumOfSamples+1);
    SonarEchosMat = zeros(2,NumOfSamples+1);
    SonarEchosMat(1,:) = NoiseLevel.*ones(1,NumOfSamples+1);
    ObsEchosVec = NoiseLevel.*ones(1,NumOfSamples+1);
    PreyEchosVec = NoiseLevel.*ones(1,NumOfSamples+1);

    
    %%% Pulse Tranmitting %%%
    for nPulseNumber = 1: NumOfTimesSonarTransmits
        CurrPulseStart = TransmittedPulsesStruct(nPulseNumber).StartPulseTime;
        CurrPulseWidth = TransmittedPulsesStruct(nPulseNumber).PulseWidth;
        CurrentFreqs = TransmittedPulsesStruct(nPulseNumber).PulseFreqsCommands;
        CurrPulseGridedFreqs = round(CurrentFreqs ./ FreqRes).* FreqRes;
        PulsePower = 10.^( TransmittedPulsesStruct(nPulseNumber).PulsePower/10);
        for kFromStartPulse= 1:CurrPulseWidth
            if CurrPulseStart+kFromStartPulse-1 <= NumOfSamples % Protecting the end time
                CurrTime = min(NumOfSamples ,CurrPulseStart+kFromStartPulse-1);
%                 try
                    SonogramMatrix(CurrPulseGridedFreqs(kFromStartPulse) , (CurrTime)) = PulsePower;
%                 catch
%                     pop=' BuildSonogram line 65' 
%                 end
                SonarEchosMat(1, CurrTime) = PulsePower;
                SonarEchosMat(2, CurrTime) = CurrentFreqs(kFromStartPulse);
            end %while CurrPulseStart+kFromStartPulse-1
        end % for kFromStartPulse= 0:CurrPulseWidth

    end % for nPulseNumber = 1: TotalNumberOfPulses
    
    %%% Echos from Obs
%     try
     calc_obs_echoes_flg = 1;
     if isfield(AllParams.SimParams,'TestMode')
         if strcmp(AllParams.SimParams.TestMode, 'swarm')
           calc_obs_echoes_flg = 0;
         end
     end % isfield(AllParams.SimParams,'TestMode')
     
     if calc_obs_echoes_flg
        [SonogramMatrix, SonarEchosMat, ObsEchosVec] = EchosSonogramAssignment(NumOfTimesObsFinds, EchosFromObsStruct, PulseWidthVec, ...
            SonarEchosMat, SonogramMatrix, ObsEchosVec, FreqRes, TransmittedPulsesStruct, AllParams );
     end % if calc_obs_echoes_flg
%     catch
%         pop = ' BuildSonogram line 80  EchosSonogramAssignment'
%     end
    
    
    %%% Echos from Prey
%     try
    [SonogramMatrix, SonarEchosMat, PreyEchosVec] = EchosSonogramAssignment(NumOfTimesPreyFinds, EchosFromPreyStruct, PulseWidthVec, ...
                SonarEchosMat, SonogramMatrix, PreyEchosVec, FreqRes, TransmittedPulsesStruct, AllParams );
%      catch
%         pop = ' BuildSonogram line 80  EchosSonogramAssignment'
%     end
    
    end % function 


   



%%
%%% EchosSonogramAssignment %%%%
function [SonogramMat, SonarEchosMat, EchosVec] = EchosSonogramAssignment(NumOfTimesFinds, EchosStruct, PulseWidthVec, ...
                SonarEchosMatrix, SonogramMatrix, EchosVector, FreqRes, TransmittedPulsesStruct, AllParams )

            
% NumOfTimesFinds = NumOfTimesObsFinds/ NumOfTimesPreyFinds
% EchosStruct = EchosFromObsStruct/ EchosFromPreyStruct
    SonarEchosMat = SonarEchosMatrix;
    SonarEchosMat1 = SonarEchosMat;
    SonogramMat = SonogramMatrix;
    EchosVec = EchosVector;
    NumOfSamples = length(EchosVec);
    SampleTime = AllParams.SimParams.SampleTime;
    NoiseLevel = 10.^(AllParams.BatSonarParams.NoiseLeveldB/10);
    %
    %% NEW CODE with 'Online' echoes
    if NumOfTimesFinds > 0 && any([EchosStruct.NumOfEchos]) % if there ara no finds - do nothing ...
        
        MaxInterVecSize= int32(AllParams.SimParams.TotalPreysNumber * NumOfTimesFinds * AllParams.BatSonarParams.PulseTimeLong/SampleTime) *2;
        FullEchosTimes = zeros(1,MaxInterVecSize);
        FullEchosPower = zeros(1,MaxInterVecSize);
        FullEchosGridedFreq = zeros(1,MaxInterVecSize);
        
        CurrIndex = 0;
        for nPulseNumber = 1: NumOfTimesFinds
            CurrPulsePower = 10.^(TransmittedPulsesStruct(EchosStruct(nPulseNumber).TransmittedPulseNum).PulsePower/10);
            CurrPulseWidth1 = EchosStruct(nPulseNumber).PulseDuration;
            CurrNumOfEchos1 = EchosStruct(nPulseNumber).NumOfEchos;
            %% new swarm 
            if CurrNumOfEchos1 > 0 % new
                CurrEchosTimes1 = [EchosStruct(nPulseNumber).EchoDetailed.EchosnTimesFull];
                CurrEchosPower = CurrPulsePower .* [EchosStruct(nPulseNumber).EchoDetailed.EchoAttenuationFull];
                CurrentFreqs1 = [EchosStruct(nPulseNumber).EchoDetailed.EchosFreqsFull];
                CurrPulseGridedFreqs1 = round(CurrentFreqs1 ./ FreqRes).* FreqRes;
                
                %%% Concatenated Vectors
                try
                    CurrIndex = CurrIndex(end)+1:(CurrIndex(end) + CurrPulseWidth1*CurrNumOfEchos1);
                    
                    %         CurrIndex = ((nPulseNumber-1)*CurrPulseWidth1*CurrNumOfEchos1+1):(nPulseNumber*CurrPulseWidth1*CurrNumOfEchos1);
                    FullEchosTimes(CurrIndex) = CurrEchosTimes1;
                    FullEchosPower(CurrIndex) = CurrEchosPower;
                    FullEchosGridedFreq(CurrIndex) = CurrPulseGridedFreqs1;
                catch
                    nPulseNumber
                    hop = ['pop nPulseNumber=', nPulseNumber]
                    %             CurrIndex
                    %             FullEchosTimes(CurrIndex)
                    %             FullEchosPower(CurrIndex)
                    %             FullEchosGridedFreq(CurrIndex)
                end % try
            end % if CurrNumOfEchos1 > 0 % new
        end % for nPulseNumber = 1: NumOfTimesFinds
        
        %%% Unique Times
        FullEchosTimes = nonzeros(FullEchosTimes)';
        UniqueEchosTimes = unique(FullEchosTimes);
        UniqueEchosTimes = min(UniqueEchosTimes,NumOfSamples); % Protecting the end samples
        PowerSum = zeros(size(UniqueEchosTimes));
        PowersMat = zeros(size(SonogramMat));
        ConcatenatedEchoesStruct(length(UniqueEchosTimes)) = ...
            struct( 'Times',[],...
            'GridedFreqs',[],...
            'MedFreq',[],...
            'Powers',[] );
        
        for kk = 1:length(UniqueEchosTimes)
            CurrTimes = [FullEchosTimes == UniqueEchosTimes(kk)];
            PowerSum(kk) = sum( FullEchosPower(CurrTimes) );
            try
                ConcatenatedEchoesStruct(kk).Times = UniqueEchosTimes(kk);
                ConcatenatedEchoesStruct(kk).GridedFreqs = [ConcatenatedEchoesStruct(kk).GridedFreqs, FullEchosGridedFreq(CurrTimes)];
                ConcatenatedEchoesStruct(kk).MedFreq = median(ConcatenatedEchoesStruct(kk).GridedFreqs );
                ConcatenatedEchoesStruct(kk).Powers = [ConcatenatedEchoesStruct(kk).Powers, FullEchosPower(CurrTimes)];
                PowersMat(ConcatenatedEchoesStruct(kk).GridedFreqs, ConcatenatedEchoesStruct(kk).Times) = ...
                    ConcatenatedEchoesStruct(kk).Powers;
            catch
                hopa ='strct Echoes'
            end% try
        end % for kk
        
        SonarEchosMat(1,[ConcatenatedEchoesStruct.Times]) =  SonarEchosMat1(1,[ConcatenatedEchoesStruct.Times]) + PowerSum;
        SonarEchosMat(2,[ConcatenatedEchoesStruct.Times]) = [ConcatenatedEchoesStruct.MedFreq];
        GridedFreqs = nonzeros([ConcatenatedEchoesStruct.GridedFreqs]);
        Times = nonzeros([ConcatenatedEchoesStruct.Times]);
        SonogramMat = SonogramMat + PowersMat;
        EchosVec(1,[ConcatenatedEchoesStruct.Times]) = EchosVec(1,[ConcatenatedEchoesStruct.Times]) + PowerSum;
        
    end %if NumOfTimesFinds > 0 % if there ara no finds - do nothing ...

end % function [SonogramMat, SonarEchosMat] = EchosSonogramAssignment(