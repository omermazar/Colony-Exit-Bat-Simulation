%% Search
% During the search the parameters are constant
Params.BatSonar.IPI_Search = 100; %msec
Params.BatSonar.TerminalFreq = 25; %kHz 
Params.BatSonar.ChirpSpan_Search = 3; % The Bandwidth of the cakk (start freq- end freq ) (kHz)
Params.BatSonar.PulseDuration_Search = 10; % msec

%% approach - 
% THe approach is defined by the distace to the taget (0.5 - 1.2 meteres)
% During the approach the bats change its paramteres linearly
Params.BatSonar.IPI_ApproachStart = 80;
Params.BatSonar.IPI_ApproachEnd = 20;
Params.BatSonar.PulseDuration_ApprachStart = 7;
Params.BatSonar.PulseDuration_ApproachEnd = 2;
Params.BatSonar.ChirpSpan_ApproachStart = 4;
Params.BatSonar.ChirpSpan_ApproachEnd = 5;

%% Buzz
% Distance is below 0.5 meters, changing parameers
% Final Buzz - distance is below 0.2 meters, Constant parameters
Params.BatSonar.TerminalFreqBuzzFinal = 23.5;
Params.BatSonar.IPI_BuzzStart = 18;  
Params.BatSonar.IPI_BuzzEnd = 10;
Params.BatSonar.IPI_BuzzFinal = 9;
Params.BatSonar.PulseDuration_BuzzStart = 2;
Params.BatSonar.PulseDuration_BuzzEnd = 1.5;
Params.BatSonar.PulseDuration_BuzzFinal = 0.75;
Params.BatSonar.ChirpSpan_Buzz = 3;
Params.BatSonar.ChirpSpan_BuzzFinal = 3;

%% General
% The Variance berween the bats in the terminal freqeuncy (kHz)
Params.BatSonar.TerminalFreqVariance = 2;
Params.BatSonar.PulsePower = 130;
Params.BatSonar.MaskingBckdSIRth = 0 ;
Params.BatSonar.MaskingFwdSIRth = 5;
Params.BatSonar.PulseDetectionTH = 40;

%% Movement
Params.BatFlight.NominalVelocity = 8.5; %m/sec
