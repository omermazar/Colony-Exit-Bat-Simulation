
function []= MyObsOffLineTest(BatDATA, BatNum)

Ts = BatDATA.AllParams.SimParams.SampleTime;

for k =BatNum:BatDATA.BAT(BatNum).NumOfTimesSonarTransmits
    t(k) = BatDATA.BAT(BatNum).TransmittedPulsesStruct(k).StartPulseTime * Ts ;
    stages(k) = string(BatDATA.BAT(BatNum).ManueverCmdStruct(k).ManueverStage);
    Type(k) = string(BatDATA.BAT(BatNum).ManueverCmdStruct(k).ManueverType);
    Power(k) = string(BatDATA.BAT(BatNum).ManueverCmdStruct(k).ManueverPower);
    DirectionCommand(k) = string(BatDATA.BAT(BatNum).ManueverCmdStruct(k).ManDirectionCommand);
    ManAccelCommand(k) = string(BatDATA.BAT(BatNum).ManueverCmdStruct(k).ManAccelCommand);

end % for

figure
plot(t, categorical(DirectionCommand),'o-')
hold on
plot(t, categorical(stages),'o-')
plot(t, categorical(Power),'o-')
plot(t, categorical(ManAccelCommand),'o-')
xlabel('time')