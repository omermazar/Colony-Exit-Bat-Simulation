function [jamPreySignal] = generateAcousicPreyJamming(preyParams, fs)

% [jamPreySignal] = generateAcousicPreyJamming(preyParams, fs)  
%                   returns the acoustic signal of the jamming moth
%                   categorized by the struct preyParams in the sampligng
%                   frequency fs (Hz).
% preyParams = AllParams.PreyParams or a struct with the followiming paramters: 
%   preyParams.
%         React2BatsFlag	    1	    % Default is 0. for Jaming moth chand to 1 [0/1]
%         DecisionPeriod	    0.01	% How long between Prey decisions? [sec]
%         JamSigDetectionTH	    30	    % the minimum signal lvl the moth can setect [dB-sPL]
%         JamBatCallsNum	    2	    % The minimum number of received bats' calls befor jamming
%         JamIPItoReact	        40      % The Maximal IPI between bats' calls to react [ms]
%         JamTxLvl	            80      % The tramsnit level of the moth's jamming signal [dB-SPL]
%         JamDelay	            0       % The time between decision and action [ms]
%         JamSwpRate	        4	    % click sweep rate [kHz/ms]
%         JamDutyCycle	        0.45    % Overall Duty-cycles of the clicks [ratio]
%         JamClickDur	        0.3	    % the duration of each click [ms]
%         JamActiveDur	        12	    % the ative cicle totalt period [ms]
%         JamTotDur	            29	    % the passive time period [ms]
%         JamISI	            5	    % Interavla between cycles [ms]
%         JamClickNum	        22	    % total number of clicks during active cycle
%         JamMainFreq	        49	    % Main frequency during the Cycle [kHz]
%         JamActiveStartFreq	82	    % The start frequency [kHz]
%         JamActiveEndFreq	    21	    % The frequency [kHz]
%         JamBW	                11	    frequency bw std in each click [kHz]
%         JamCyclesNum	        2	    active-active- passive-passive
%         Jam2ndStart	        4	    4the second active start time [ms]
%         JamRepetitions        5       Number of Consecutive repetitions when trigerred

%% Genearl Parameters 
periodNsamples    = round(preyParams.JamTotDur * 1e-3 * fs);
if preyParams.JamCyclesNum >1
    periodNsamples = periodNsamples + preyParams.Jam2ndStart * 1e-3 * fs;
end
jamNsamples = periodNsamples;
% repeaptitions
% if preyParams.JamRepetitions > 1
%     jamNsamples = periodNsamples*(preyParams.JamRepetitions - 1);
% else
%     jamNsamples = periodNsamples;
% end
jamPreySignal  = zeros(1, jamNsamples);   % output

activeNsamples  = round(preyParams.JamActiveDur * 1e-3 * fs);
passiveNsamples = activeNsamples;

% Clicks
numberOfClicks = preyParams.JamClickNum;
clickNSamples   = round(preyParams.JamClickDur * 1e-3 * fs);
silenceNsamples = round((activeNsamples -  clickNSamples * numberOfClicks ) ./ (numberOfClicks -1) );

clickBw   = preyParams.JamBW * 1e3;
freqHigh = preyParams.JamActiveStartFreq * 1e3 - clickBw; %Hz
freqLow   = preyParams.JamActiveEndFreq* 1e3 + clickBw; %Hz
freqStep  = (freqHigh - freqLow) ./ (numberOfClicks-1);

%% The Active Cycle (#1)
startIx = 1;
activeIx = startIx:startIx + activeNsamples-1;

[activeSignal] = ClickTrain( freqHigh, freqStep, activeNsamples, clickNSamples, clickBw, silenceNsamples, numberOfClicks, fs);
jamPreySignal(activeIx) = jamPreySignal(activeIx) + activeSignal;

%%  The Passive Cycle #1
startPassiveIx = activeNsamples + preyParams.JamISI * 1e-3 * fs ;
passiveIx = startPassiveIx:startPassiveIx + activeNsamples -1;

passiveSignal = fliplr(activeSignal);
jamPreySignal(passiveIx) = jamPreySignal(passiveIx) + passiveSignal;

% [jamPreySignal] = ClickTrain( freqLow, -freqStep, clickNSamples, clickBw, silenceNsamples, numberOfClicks, fs);

%% The Active Cycle #2
startIx = 1 + preyParams.Jam2ndStart * 1e-3 *fs ;
activeIx = startIx:startIx + activeNsamples -1;

jamPreySignal(activeIx) = jamPreySignal(activeIx) + activeSignal;
% [jamPreySignal] = ClickTrain(jamPreySignal, startIx, freqHigh, freqStep, clickNSamples, clickBw, silenceNsamples, numberOfClicks, fs);

%%  The Passive Cycle #2
startPassiveIx = activeNsamples + preyParams.Jam2ndStart * 1e-3 *fs  + preyParams.JamISI * 1e-3 * fs ;
passiveIx = startPassiveIx:startPassiveIx + activeNsamples - 1;
jamPreySignal(passiveIx) = jamPreySignal(passiveIx) + passiveSignal;

% [jamPreySignal] = ClickTrain(jamPreySignal, startPassiveIx, freqLow, -freqStep, clickNSamples, clickBw, silenceNsamples, numberOfClicks, fs);

%% repeatitions
if preyParams.JamRepetitions >1
    % padding
    jamPreySignal((end+1):(end+preyParams.Jam2ndStart * 1e-3 * fs) ) = 0;
    startRepIx = periodNsamples + preyParams.Jam2ndStart * 1e-3 * fs ;
    jamPreySignal = repmat(jamPreySignal, 1, preyParams.JamRepetitions );
%     for k = 1:(preyParams.JamRepetitions - 1)
%         periodIx = startRepIx:startRepIx + periodNsamples - 1;
%         jamPreySignal(periodIx) = jamPreySignal(1:periodNsamples);
%         startRepIx =  startRepIx + preyParams.Jam2ndStart * 1e-3 * fs ;
%     end % for k
end % if preyParams.JamRepetitions >1




%% Ouprput Lvl
sigSTD = std(jamPreySignal);
jamPreySignal = jamPreySignal ./ sigSTD * 10.^(preyParams.JamTxLvl / 20);



end %  main function

%% Function Click Train generation
function [sig_out] = ClickTrain(freqStart, freqStep, cycleNsamples, clickNSamples, clickBw, silenceNsamples, numberOfClicks, fs )
minFreq = 20e3;
maxFreq = 80e3;

% sig_out = sig_in; 
sig_out = zeros(1,cycleNsamples);
startIx = 1;

for n = 0:numberOfClicks-1
    currFreq = freqStart - n * freqStep; 
    fhigh = currFreq + clickBw/2;
    flow  = currFreq - clickBw/2;
    t = 0: 1/fs : clickNSamples/fs;
%     x1 =  chirp(t, fhigh, t(end), flow,'logarithmic');
    
    rng(n)   
    whN = randn(1, clickNSamples+1);
    d = designfilt('bandpassiir', 'FilterOrder',12, ...
        'HalfPowerFrequency1', flow,'HalfPowerFrequency2', fhigh, ...
        'SampleRate',fs);
    x1 = filter(d,whN);
    
    currIx = startIx + (0 : clickNSamples);
    sig_out(currIx) = sig_out(currIx) + x1;
    
    % add lower harmonic
    if currFreq/2 > minFreq
        flow  = max(currFreq/2 - clickBw/3, minFreq);
        fhigh = min(currFreq/2 + clickBw/3, maxFreq);

        d = designfilt('bandpassiir', 'FilterOrder',12, ...
            'HalfPowerFrequency1', flow,'HalfPowerFrequency2', fhigh, ...
            'SampleRate',fs);
        x1 = 0.5 * filter(d,whN);
        sig_out(currIx) = sig_out(currIx) + x1;
    end
    
    % add upper harmoinc
    % add lower harmonic
    if currFreq*2 < maxFreq
        flow  = max(currFreq*2 - clickBw/3, minFreq);
        fhigh = min(currFreq*2 + clickBw/3, maxFreq);

        d = designfilt('bandpassiir', 'FilterOrder',12, ...
            'HalfPowerFrequency1', flow,'HalfPowerFrequency2', fhigh, ...
            'SampleRate',fs);
        x1 = 0.5 * filter(d,whN);
        sig_out(currIx) = sig_out(currIx) + x1;
    end
    startIx = startIx + clickNSamples + silenceNsamples ;
end % for n

end % function