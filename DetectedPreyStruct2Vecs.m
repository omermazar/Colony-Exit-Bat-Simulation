function [VecDistances, VecAngles ] = DetectedPreyStruct2Vecs(Prey) 

NumOfPreyDetected = Prey.NumOfPreyDetected;
DetectedPreys = Prey.DetectedPreys;
VecDistances = zeros(1, NumOfPreyDetected);
VecAngles = zeros(1, NumOfPreyDetected);
for n = 1:NumOfPreyDetected
    VecDistances(n) = Prey.Prey(DetectedPreys(n)).Distance;
    VecAngles(n) = Prey.Prey(DetectedPreys(n)).Alpha;
end 
