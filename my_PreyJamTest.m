function [] = my_PreyJamTest(BatDATA, preyNum, batNum) 

% plot the jamming times and places of the moths
fig = figure;
axXY = subplot(2,1,1); hold(axXY,'on')
axTime = subplot(2,1,2); hold(axTime,'on')
fig.Position = [400 150 560 600 ];

prey = BatDATA.PREY(preyNum);
bat = BatDATA.BAT(batNum);

%% the XY PLOt
% prey position
plot(axXY, prey.xPrey, prey.yPrey, '.b', 'DisplayName','Prey')

% the PreyJammingLocations
ixJam = round(prey.JammingnTimes);
plot(axXY, prey.xPrey(ixJam), prey.yPrey(ixJam), '*r', 'DisplayName','Jamming')

% the Bat Postion
plot(axXY, bat.xBatPos, bat.yBatPos, '.', 'DisplayName','Bat');
title(axXY, ['Jamming Moths: Prey ', num2str(preyNum), ', Bat ', num2str(batNum)])
lXY = legend(axXY);
axXY.Position = [0.13 0.43 0.775 0.5];
lXY.Location = 'bestoutside';

%% the time plot
[sortedTimes, ixSorted ] = sort(prey.RxnTimes);
% plot(axTime, prey.RxnTimes, 10*log10(prey.RxLvls), '.k', 'DisplayName', 'Rxlvl all')
plot(axTime, sortedTimes, 10*log10(prey.RxLvls(ixSorted)), '.-k', 'DisplayName', 'Rxlvl all')
plot(axTime, prey.RxnTimes(prey.RxBatNumTx == batNum), 10*log10(prey.RxLvls(prey.RxBatNumTx == batNum)), ...
    'ob', 'DisplayName', ['BAT', num2str(batNum)])
rxTH = BatDATA.AllParams.PreyFlightParams.JamSigDetectionTH;
plot(axTime, axTime.XLim, rxTH*[1,1], 'r--', 'DisplayName', 'DetectionTH')
plot(axTime, prey.JammingnTimes, rxTH*ones(size(prey.JammingnTimes)), 'r*', 'DisplayName', 'Jamming')
ylabel(axTime, 'dB')
grid(axTime, 'on')
lTime = legend(axTime);
xlabel(axTime, 'samples')
axTime.Position = [0.13 0.11 0.775 0.28];
lTime.Location = 'bestoutside';




