function [DFErr] = calcDFError(SNR, inAngles)

% Anoise = d/lambda, d= head widfth~2cm, lamda@40kHz = 8mm
Anoise = 0.5*pi/180;% *10; % 1deg error at SNR =10db
% DFNoiseErrVec =  min (Anoise / sqrt(NumberOfMeasures) ./ SNR , pi/2 ) ; % the max error is pi/2
DFNoiseErrVec =  min (Anoise ./ SNR , pi/2 ) ; % the max error is pi/2

% the resolution of the measurement as function of teta

% blur = 0.5deg @angle= 0deg, ~10deg~angle= 90deg
a = 1*pi/180; 
b = 0.2*pi/180; % 1.5 old
DFBlurrVec = b + a.*abs(sin(inAngles));

DFErrSTD = sqrt(DFNoiseErrVec.^2 + DFBlurrVec.^2);
DFErr = DFErrSTD .* randn(size(inAngles));