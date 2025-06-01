figure; hold on; grid on

nCalls = 30;
for n = 10:nCalls
    %% refence
    t1 = BatDATA_Feb08.BAT(1).TransmittedPulsesStruct(n).StartPulseTime;
    plot(BatDATA_Feb08.BAT(1).ObsFindsStruct(n).xFindsEst, BatDATA_Feb08.BAT(1).ObsFindsStruct(n).yFindsEst, 'og')
    quiver(BatDATA_Feb08.BAT(1).xBati(t1), BatDATA_Feb08.BAT(1).yBati(t1), ...
        50*cos(BatDATA_Feb08.BAT(1).Teta(t1) + BatDATA_Feb08.BAT(1).ManueverCmdStruct(n).Angle2HuntedPrey), ...
        50*sin(BatDATA_Feb08.BAT(1).Teta(t1) + BatDATA_Feb08.BAT(1).ManueverCmdStruct(n).Angle2HuntedPrey), 'o-g')
    
    %% OJT
    t2 = BatDATA_OJT.BAT(1).TransmittedPulsesStruct(n).StartPulseTime;
    plot(BatDATA_OJT.BAT(1).ObsFindsStruct(n).xFindsEst, BatDATA_OJT.BAT(1).ObsFindsStruct(n).yFindsEst, '.k')
    quiver(BatDATA_OJT.BAT(1).xBati(t2), BatDATA_OJT.BAT(1).yBati(t2), ...
        50*cos(BatDATA_OJT.BAT(1).Teta(t2) + BatDATA_OJT.BAT(1).ManueverCmdStruct(n).Angle2HuntedPrey), ...
        50*sin(BatDATA_OJT.BAT(1).Teta(t2) + BatDATA_OJT.BAT(1).ManueverCmdStruct(n).Angle2HuntedPrey), 'o-k')

    %% DeBUg
    t3 = BatDATA01.BAT(1).TransmittedPulsesStruct(n).StartPulseTime;
    plot(BatDATA01.BAT(1).ObsFindsStruct(n).xFindsEst, BatDATA01.BAT(1).ObsFindsStruct(n).yFindsEst, '.b')
    quiver(BatDATA01.BAT(1).xBati(t3), BatDATA01.BAT(1).yBati(t3), ...
        50*cos(BatDATA01.BAT(1).Teta(t3) + BatDATA01.BAT(1).ManueverCmdStruct(n).Angle2HuntedPrey), ...
        50*sin(BatDATA01.BAT(1).Teta(t3) + BatDATA01.BAT(1).ManueverCmdStruct(n).Angle2HuntedPrey), 'o-b')
end

figure; subplot(2,1,1); hold on
plot(categorical(string({BatDATA_OJT.BAT(1).ManueverCmdStruct(1:nCalls).ManueverType})), 'o-k')
plot(categorical(string({BatDATA_Feb08.BAT(1).ManueverCmdStruct(1:nCalls).ManueverType})), 'o-g')
plot(categorical(string({BatDATA01.BAT(1).ManueverCmdStruct(1:nCalls).ManueverType})), 'o-b')
title('ManType')

subplot(2,1,2); hold on
plot([BatDATA_OJT.BAT(1).ManueverCmdStruct(1:nCalls).Angle2HuntedPrey], 'o-k');
plot([BatDATA_Feb08.BAT(1).ManueverCmdStruct(1:nCalls).Angle2HuntedPrey], 'o-g');
plot([BatDATA01.BAT(1).ManueverCmdStruct(1:nCalls).Angle2HuntedPrey], 'o-b');
title('Directon Command')

% subplot(3,1,3); hold on
% plot([BatDATA_OJT.BAT(1).Teta(1:t2)], 'o-k');
% plot([BatDATA_Feb08.BAT(1).Teta(1:t1)], 'o-g');
% title('FlightDirecton')
