function [InterSonogram] = BuildInterferenceSonogramFromOnLine(BAT , AllParams, varargin)

%[InterSonogram] = BuildInterferenceSonogram(BAT , AllParams)
% will build the sonogram of the interfence from other bats to current bat
%
% [InterSonogram] = BuildInterferenceSonogram(BAT , AllParams, BatSonogram)
% Will add the interference to cuurent BAT sonogram


% Consts and Inputs
NumOfBats = AllParams.SimParams.TotalBatsNumber;
SimulationTime = AllParams.SimParams.SimulationTime;
SampleTime = AllParams.SimParams.SampleTime;
NumOfSamples = SimulationTime/SampleTime;
MaxFreq = AllParams.BatSonarParams.CenterFreq + ...
    AllParams.BatSonarParams.ChirpSpan*AllParams.BatSonarParams.ChirpFlag*0.5;
FreqRes = 1;
FreqGrid = unique(round([1:MaxFreq+10]/FreqRes))*FreqRes;

NoiseLevel = 10.^(AllParams.BatSonarParams.NoiseLeveldB/10);

if nargin == 3
    InterSonogram = varargin{1};
else
    InterSonogram = NoiseLevel.*ones(length(FreqGrid) ,NumOfSamples+1);
end

% START
InterferenceIndex = find(abs(BAT.InterferenceVec) ~= NoiseLevel);
nTimes = min([BAT.InterferenceFullStruct(InterferenceIndex).Times] , NumOfSamples+1);

for nT = nTimes
    CurrFreqs = BAT.InterferenceFullStruct(nT).Freqs;
    CurrGridedFreqs = abs(round(CurrFreqs ./ FreqRes .* FreqRes(1)));
    CurrGridedFreqs = min(CurrGridedFreqs, size(InterSonogram,1));
    CurrPowers = BAT.InterferenceFullStruct(nT).Power(1);
    try
        InterSonogram(CurrGridedFreqs, nT) = InterSonogram(CurrGridedFreqs, nT) + CurrPowers';
    catch
        pop = ['fuck', num2str(nT)]
    end % try
end % for nT




% for kBAT = 1:NumOfBats
%    kBatFreqs = BAT.BatInterference(kBAT).InterFreqs;
%    kBatGridedFreqs = abs(unique(round(kBatFreqs ./ FreqRes .* FreqRes) ));
%    TimesOfInterference = BAT.BatInterference(kBAT).RecivedInterPulseTimes;
%    InterPower = BAT.BatInterference(kBAT).InterPower;
%    for nn= 1:length(TimesOfInterference)
%        CurrFreq = abs(round(kBatFreqs(nn) ./ FreqRes .* FreqRes));
%        try
%            InterSonogram(CurrFreq , TimesOfInterference(nn)) = ...
%                InterSonogram(CurrFreq , TimesOfInterference(nn)) + InterPower(nn);
%        catch
%            pop0= 'booo  function BuildInterferenceSonogram  line 40'
%        end %try
%    end %for nn= 1:length(TimesOfInterference)
% end  %for kBAT = 1:NumOfBats
    


