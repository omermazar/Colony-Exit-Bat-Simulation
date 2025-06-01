function[ MissDetectedByPD ]= IsPreyMissDetected(DetectedPreys, ProbToDetect)

% IsPreyMissDetected(IsPreyMissDetected(DetectedPreys, ProbToDetect) 
%
% Returns the pulses that are miised by Probalbilty of detection

NumOfDetections = length(DetectedPreys);
ProbVec = rand(1,NumOfDetections);
MissDetectedByPD = DetectedPreys(find(ProbVec > ProbToDetect));

