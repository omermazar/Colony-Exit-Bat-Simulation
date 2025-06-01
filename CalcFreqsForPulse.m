function [FreqsVec, ChirpMinMaxFreqs] = CalcFreqsForPulse(PulseWidth, MinFreq, ChirpSpan, ChirpShape)

% this function will calculate the frequncies of the pulse
%  FreqsVec = CalcFreqsForPulsength (PulseWidth, CenterFreq, ChirpSpan, ChirpType)
%
% Output - 
%   FreqsVec - the vector of the freqs length of PulseDuration
%   (round(PulseWidth)
%   %  ChirpMinMaxFreqs - the min Freq and the max freq of the chirp (not
%  rounded)
% Input - 
%   PulseWidth - the length of the pulse transmitted
%   CenterFreq, ChirpSpan - parmetertes of the pulse
%   ChirpType - 'Foraging' / 'RegularHunt' /  'Buzz' / 'RegularManuver' /
%   'CrushAvoidance' 
%   ChirpShape = 'Linear' or 'Geometric' or 'Logarithmic''


PulseDuration = max(1, round(PulseWidth));

if PulseDuration >2 
    MaxFreq = MinFreq + ChirpSpan;
else %if PulseDuration >2 ,  if the is short, cut the chirp down
    MaxFreq = MinFreq + 0.5*ChirpSpan;
end % if PulseDuration >2 

ChirpMinMaxFreqs = [MinFreq, MaxFreq];

% FreqsVec = zeros(1,PulseWidth);
if PulseDuration == 1
    FreqsVec = MinFreq;
else % if PulseDuration == 1
    Den = max( PulseDuration-1 , 1);
    switch ChirpShape
        case 'Linear'
            D = (MaxFreq - MinFreq)/Den;
            FreqsVec = [MaxFreq : -D : MinFreq];
        case 'Geometric' 
            Q = (MinFreq/MaxFreq)^(1/Den);
            FreqsVec = MaxFreq * Q.^(0:Den);
        case 'Logarithmic'
              A = MaxFreq;
              b = log10(MinFreq/MaxFreq)./PulseDuration;
              FreqsVec = A*10.^(b*(1:PulseDuration));
%             A = (MaxFreq - MinFreq)/log10(PulseDuration+.01); % avoiding(log(1)
%             FreqsVec = MaxFreq - A.*log10(1:PulseDuration);
%               t=0:PulseDuration;
%               Freq2 =  chirp(t, MaxFreq, max(t), MinFreq,'logarithmic');
    end % switch ChirpShape
end % if PulseDuration == 1


