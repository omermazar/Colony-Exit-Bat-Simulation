function [T] = getEchoDataforTable(EchoesStruct, TxBat, echoType, txPower, maxDelay, AllParams)

%%%% New June2024
% this function returne the relevant data for the full echoes table


nEchoes = EchoesStruct.NumOfEchos;
T.time = (EchoesStruct.TransmittedPulseTime + EchoesStruct.EchosTimes)' * AllParams.SimParams.SampleTime;
T.callNum = EchoesStruct.TransmittedPulseNum*ones(nEchoes,1);
T.distance = EchoesStruct.EchosmDistances';
T.directionOfArrival = EchoesStruct.EchosAngles';
T.receivedLevel = EchoesStruct.EchosAttenuation' .* txPower;
T.TxBat = repmat(TxBat, nEchoes,1);
T.echoType = repmat(echoType, nEchoes,1);
T.targetIndex = EchoesStruct.TargetIndex';
T = struct2table(T);
% Only echoes above TH level
% ix = T.receivedLevel >  10.^(AllParams.BatSonarParams.PulseDetectionTH/10);
% Consider echoes from current Call Only (no Ambiguity)
ix = T.time <= EchoesStruct.TransmittedPulseTime.*AllParams.SimParams.SampleTime + maxDelay;
T = T(ix,:);

