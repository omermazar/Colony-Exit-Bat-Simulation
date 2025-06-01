function [InterSonogram] = BuildInterferenceSonogram(BAT , AllParams, varargin)

%[InterSonogram] = BuildInterferenceSonogram(BAT , AllParams)
% will build the sonogram of the interfence from other bats to current bat
%
% [InterSonogram] = BuildInterferenceSonogram(BAT , AllParams, BatSonogram)
% Will add the interference to cuurent BAT sonogram


% Consts an Inputs
NumOfBats = AllParams.SimParams.TotalBatsNumber;
SimulationTime = AllParams.SimParams.SimulationTime;
SampleTime = AllParams.SimParams.SampleTime;
NumOfSamples = SimulationTime/SampleTime;
MaxFreq = AllParams.BatSonarParams.CenterFreq + ...
    AllParams.BatSonarParams.ChirpSpan*AllParams.BatSonarParams.ChirpFlag*0.5;
% % % PulsePower = AllParams.BatSonarParams.PulsePower;
FreqRes = 1;
FreqGrid = unique(round([0:MaxFreq]/FreqRes))*FreqRes;

% NoiseLevel = 0.1;
NoiseLevel = 10.^(AllParams.BatSonarParams.NoiseLeveldB/10);

if nargin == 3
    InterSonogram = varargin{1};
else
    InterSonogram = NoiseLevel.*ones(length(FreqGrid) ,NumOfSamples+1);
end

for kBAT = 1:NumOfBats
   kBatFreqs = BAT.BatInterference(kBAT).InterFreqs;
   kBatGridedFreqs = abs(unique(round(kBatFreqs ./ FreqRes .* FreqRes) ));
   TimesOfInterference = BAT.BatInterference(kBAT).RecivedInterPulseTimes;
   InterPower = BAT.BatInterference(kBAT).InterPower;
   for nn= 1:length(TimesOfInterference)
       CurrFreq = abs(round(kBatFreqs(nn) ./ FreqRes .* FreqRes));
%        try
           InterSonogram(CurrFreq , TimesOfInterference(nn)) = ...
               InterSonogram(CurrFreq , TimesOfInterference(nn)) + InterPower(nn);
%        catch
%            pop0= 'booo  function BuildInterferenceSonogram  line 40'
%        end %try
   end %for nn= 1:length(TimesOfInterference)
end  %for kBAT = 1:NumOfBats
    


