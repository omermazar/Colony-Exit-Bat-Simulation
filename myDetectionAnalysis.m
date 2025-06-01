function [stats] = myDetectionAnalysis(BatDATA, detDistIn, obsDet)
% calculate the probabilty of detecions and histogram of distances 

%% Init
detAllDist  = [];
distAll     = [];
anglesAll   = [];
detAllAngle = [];
detAllRxLvl = [];

ownCallsAllLvl       = [];
preyAllEchoesLvl     = [];
conpsAllEchoesLvl    = [];
obsAllEchoesLvl      = [];
conspsObsEchosLvl    = [];
conspsConspsEchosLvl = [];
conpsCallsAllLvl     = [];
totalInterlLvl       = [];

if nargin < 3
    obsDet = 'Conspecific';
end


%% data



%%%% Grab data of all bats
for k=1:BatDATA.AllParams.SimParams.TotalBatsNumber
    % the detections
    switch obsDet
        case 'Conspecific'
            detAllDist = [detAllDist, [BatDATA.BAT(k).Consps_FindsStruct.Distances]];
            detAllAngle = [detAllAngle, [BatDATA.BAT(k).Consps_FindsStruct.Angles]];
            % the distances
            distAll   = [distAll, BatDATA.BAT(k).OtherBatsPolarCoordinates.Distances];
            anglesAll = [anglesAll, BatDATA.BAT(k).OtherBatsPolarCoordinates.TargetAngle];
            % SignalLevels
            detAllRxLvl = [detAllRxLvl, [BatDATA.BAT(k).Consps_FindsStruct.RxPowerOfDetectedPreys]];

        case 'Prey'
            detAllDist = [detAllDist, [BatDATA.BAT(k).PreyFindsStruct.Dist2DetectedPrey]];
            detAllAngle = [detAllAngle, [BatDATA.BAT(k).PreyFindsStruct.Angle2DetectedPrey]];
            % the distances
            currDist  = [BatDATA.BAT(k).Vector2Preys.Dist2Prey];
            currAngle = [BatDATA.BAT(k).Vector2Preys.Angle2Prey];
            % find the times of calls by the derivative
            ixCalls = diff(currDist) ~= 0;
            distAll   = [distAll, currDist(ixCalls)];
            anglesAll = [anglesAll, currAngle(ixCalls)];
            % SignalLevels
            detAllRxLvl = [detAllRxLvl, [BatDATA.BAT(k).PreyFindsStruct.RxPowerOfDetectedPreys]];

        case 'Obs'
            detAllDist = [detAllDist, [BatDATA.BAT(k).ObsFindsStruct.Distances ]];
            detAllAngle = [detAllAngle, [BatDATA.BAT(k).ObsFindsStruct.Angles]];
            % the distances
            distAll   = [distAll, BatDATA.BAT(k).ObsInBeamStruct.Distances];
            anglesAll = [anglesAll, BatDATA.BAT(k).ObsInBeamStruct.TargetAngle];
            % SignalLevels
            detAllRxLvl = [detAllRxLvl, [BatDATA.BAT(k).ObsFindsStruct.RxPowerOfDetectedPreys]];

    end % switch
    
    % Omit Inf 
    distAll(distAll == inf) = nan;
    anglesAll(distAll == inf) = nan;

    %%%%% Acoustics
    ownCallsAllLvl = [ownCallsAllLvl, [BatDATA.BAT(k).TransmittedPulsesStruct.PulsePower] ] ;

    if BatDATA.AllParams.SimParams.AcousticsCalcFlag

          % All Prey
          ix = 20*log10(abs(BatDATA.BAT(k).AcousticSig_PreyEchoes)) > 1;
          preyAllEchoesLvl = [preyAllEchoesLvl, 20*log10(abs(BatDATA.BAT(k).AcousticSig_PreyEchoes(ix)))] ;
          % All Conspsefics Echoes
          ix = 20*log10(abs(BatDATA.BAT(k).AcousticSig_CospsEchoes)) > 1;
          conpsAllEchoesLvl = [conpsAllEchoesLvl, 20*log10(abs(BatDATA.BAT(k).AcousticSig_CospsEchoes(ix)))] ;
          % All Obs
          ix = 20*log10(abs(BatDATA.BAT(k).AcousticSig_Clutter)) > 1;
          obsAllEchoesLvl = [obsAllEchoesLvl, 20*log10(abs(BatDATA.BAT(k).AcousticSig_Clutter(ix)))] ;
          % Obs from consps
          ix = 20*log10(abs(BatDATA.BAT(k).AcousticSig_ConspsObsEchoes)) > 1;
          conspsObsEchosLvl = [conspsObsEchosLvl, 20*log10(abs(BatDATA.BAT(k).AcousticSig_ConspsObsEchoes(ix)))] ;
          % Consps from consps
          ix = 20*log10(abs(BatDATA.BAT(k).AcousticSig_ConspsConspsEchoes)) > 1;
          conspsConspsEchosLvl = [conspsConspsEchosLvl, 20*log10(abs(BatDATA.BAT(k).AcousticSig_ConspsConspsEchoes(ix)))] ;

          % Masking Signal - Consps calls
          ix = 10*log10(abs(BatDATA.BAT(k).AcousticSig_ConspsCalls)) > 1;
          conpsCallsAllLvl = [conpsCallsAllLvl, 20*log10(abs(BatDATA.BAT(k).AcousticSig_ConspsCalls(ix)))] ;

          % Total Masking Signal
          ix = 10*log10(abs(BatDATA.BAT(k).AcousticSig_TotalInterference)) > 1;
          totalInterlLvl = [totalInterlLvl, 20*log10(abs(BatDATA.BAT(k).AcousticSig_TotalInterference(ix)))] ;
   
    else %  if BatDATA.AllParams.SimParams.AcousticsCalcFla
       % All Prey
          ix = 10*log10(abs(BatDATA.BAT(k).PreyEchosVec)) > 1;
          preyAllEchoesLvl = [preyAllEchoesLvl, 10*log10(abs(BatDATA.BAT(k).PreyEchosVec(ix)))] ;
          % All Conspsefics Echoes
          ix = 10*log10(abs(BatDATA.BAT(k).Consps_EchosVec)) > 5;
          conpsAllEchoesLvl = [conpsAllEchoesLvl, 10*log10(abs(BatDATA.BAT(k).Consps_EchosVec(ix)))] ;
          % All Obs
          ix = 10*log10(abs(BatDATA.BAT(k).ObsEchosVec)) > 1;
          obsAllEchoesLvl = [obsAllEchoesLvl, 10*log10(abs(BatDATA.BAT(k).ObsEchosVec(ix)))] ;
          % Obs from consps
%           ix =10*log10(abs(BatDATA.BAT(k).PreyEchosVec)) > 1;
%           conspsObsEchosLvl = [conspsObsEchosLvl, 10*log10(abs(BatDATA.BAT(k).PreyEchosVec(ix)))] ;
          % Masking Signal - Consps calls
          ix = 10*log10(abs(BatDATA.BAT(k).AllInterPulses)) > 1;
          conpsCallsAllLvl = [conpsCallsAllLvl, 10*log10(abs(BatDATA.BAT(k).AllInterPulses(ix)))] ;
          % Total Masking Signal
          ix = 10*log10(abs(BatDATA.BAT(k).AllInterPulses)) > 1;
          totalInterlLvl = [totalInterlLvl, 10*log10(abs(BatDATA.BAT(k).AllInterPulses(ix)))] ;
    end %  if BatDATA.AllParams.SimParams.AcousticsCalcFla
end % for
% go back to meters
detAllDist = detAllDist*0.01; %meters
if ~strcmp(obsDet, 'Prey')
    distAll    = distAll*0.01; %meters
end % if
% detection probabilty
if nargin < 2
    detDist  = 3;
else
    detDist = detDistIn;
end
detAngle = pi/3;

ixDetDetIn = detAllDist < detDist & abs(detAllAngle)<detAngle;
ixDetAllIn = distAll < detDist & abs(anglesAll)<detAngle;
% detectProb = numel(detAllDist) ./ sum(ixDetAll);
detectProbIn = sum(ixDetDetIn) ./ sum(ixDetAllIn);

% detection Probabity as a function of target distance
detDist = [0.2:0.1:3];
detectProb = nan(size(detDist));
totDetAll  = nan(size(detDist));

for kk = 1:numel(detDist)
    ixDetDet = detAllDist < detDist(kk) & abs(detAllAngle)<detAngle;
    ixDetAll = distAll < detDist(kk) & abs(anglesAll)<detAngle;
    % detectProb = numel(detAllDist) ./ sum(ixDetAll);
    totDetAll(kk) = sum(ixDetAll);
    if min([sum(ixDetDet), sum(ixDetAll)]) > 3 & sum(ixDetAll) > 0
        detectProb(kk) = sum(ixDetDet) ./ sum(ixDetAll);
    end
end

%% Output
stats.totalDetections    = numel(detAllDist);
stats.beamRange = detDist;
stats.beamAngle = detAngle;
stats.inBeamDetections   = sum(ixDetDet);
stats.inBeamTotal        = sum(ixDetAll);
stats.detectProb         = detectProb;
[N1,edges1] = histcounts(detAllDist,20);
[N2,edges2] = histcounts(distAll,20);
stats.histDetectedDist.binCounts   = N1;
stats.histDetectedDist.egdes       = edges1;
stats.histAllDist.ninCounts        = N2;
stats.histAllDist.egdes            = edges2;

stats.meanDetectedDist   = mean(detAllDist);
stats.stdDetectedDist    = std(detAllDist);
stats.medianDetectedDist = median(detAllDist);

stats.meanDetectedRxLvl   = mean(detAllRxLvl);
stats.stdDetectedRxLvl    = std(detAllRxLvl);
stats.medianDetectedRxLvl = median(detAllRxLvl);

%% Figures
nBats = BatDATA.AllParams.SimParams.TotalBatsNumber;
nPrey = BatDATA.AllParams.SimParams.TotalPreysNumber;

fig = figure;

fig.Position([2,4]) = [100, 630]; % Chanעe the height og the figure
fig.Position([3]) = [720]; % Change the width og the figure
ax1 = subplot(3,2,1); hold on
% histogram(distAll, 20, 'FaceColor', 'b', 'Normalization','probability');
% histogram(detAllDist, 20, 'FaceColor', 'g', 'Normalization','probability');
Y = {distAll, detAllDist};
fc = [0 0 1 ; 0 1 0]; % 'b';'g'
txt = {'all distances', 'detected distances' }
violin(Y, 'facecolor', fc, 'edgecolor', 'b', 'mc', 'k', 'medc', 'r--');
ax1.FontSize = 9;
ax1.XTick = 1:size(Y,2);
ax1.XTickLabel = txt;
% title
title('Distriburion of Distances ');

grid on
ylabel('Target Distance (m)')
legend('off');
% hist1 = histogram(detAllDist,20);
% title(['Distriburion of Distances ', num2str(nBats),'Bats, ', num2str(nPrey), 'Prey-Items'])

x1 = 1.7; %max(xlim)*0.5;
y1 = max(ylim)*0.8;
y2 = max(ylim)*0.35;
y3 = max(ylim)*0.7;
y4 = max(ylim)*0.25;
t1 = text(x1, y3, ['total: ', num2str(numel(detAllDist))],   'FontSize', 7);
t2 = text(x1, y1, ['mean dist: ', num2str(stats.meanDetectedDist,2), '+-',  num2str(stats.stdDetectedDist,2)],  'FontSize', 7);
% text(x1, y3, ['Beam Range: ', num2str(detDistIn)] )    
% text(x1, y4, ['detection probabity: ', num2str(detectProbIn,2)] )    

% RxPower
ax11 = subplot(3,2,2); hold on; grid on
Y = {detAllRxLvl};
txt = {'Detected'};
fc = [0 1 0]; % 'g'
% if BatDATA.AllParams.SimParams.AcousticsCalcFlag
    switch obsDet
        case 'Prey'
%             histogram(ax22, preyAllEchoesLvl, 100, 'FaceColor', 'b', 'Normalization','probability', 'DisplayName', 'All');
%             histogram(ax11, preyAllEchoesLvl, 100, 'FaceColor', 'b', 'Normalization','probability', 'DisplayName', 'Prey');
%             histogram(ax11, conpsCallsAllLvl, 100, 'FaceColor', 'r', 'Normalization','probability', 'DisplayName', 'ConspsCalls');
%             histogram(ax11, obsAllEchoesLvl, 100, 'FaceColor', 'y', 'Normalization','probability', 'DisplayName', 'Clutter');
            Y{end+1}  = preyAllEchoesLvl;
            txt{end+1} = 'AllPrey';
            fc(end+1,:) = [0 0 1]; % 'b'
            if ~isempty(conpsCallsAllLvl)
                Y{end+1}   = conpsCallsAllLvl;
                txt{end+1} = 'ConspsCalls';
                fc(end+1,:) = [1 0 0]; %'r';
            end
             if ~isempty(obsAllEchoesLvl)
                Y{end+1}   = obsAllEchoesLvl;
                txt{end+1} = 'Clutter';
                fc(end+1,:) = [1 1 0]; %'y'
             end
             if ~isempty(conspsConspsEchosLvl)
                Y{end+1}   = conspsConspsEchosLvl;
                txt{end+1} = 'ConspsByConsps';
                fc(end+1,:) = [0.5 0.5 0.5]; %'grey'
             end
             
        case 'Obs'
%             histogram(ax22, obsAllEchoesLvl, 100, 'FaceColor', 'b', 'Normalization','probability', 'DisplayName', 'All');
%             histogram(ax11, conpsCallsAllLvl, 100, 'FaceColor', 'r', 'Normalization','probability', 'DisplayName', 'ConspsCalls');
%             histogram(ax11, conspsObsEchosLvl, 100, 'FaceColor', 'y', 'Normalization','probability', 'DisplayName', 'ConspsObsEchoes');
%             histogram(ax11, conpsAllEchoesLvl, 100, 'FaceColor', 'm', 'Normalization','probability', 'DisplayName', 'ConspsEchoes');
            Y{end+1}  = obsAllEchoesLvl;
            txt{end+1} = 'AllObs';
            fc(end+1,:) = [0 0 1]; % 'b'
            if ~isempty(conpsCallsAllLvl)
                Y{end+1}   = conpsCallsAllLvl;
                txt{end+1} = 'ConspsCalls';
                fc(end+1,:) = [1 0 0]; %'r';
            end
             if ~isempty(conspsObsEchosLvl)
                Y{end+1}   = conspsObsEchosLvl;
                txt{end+1} = 'ConspsObsEchoes';
                fc(end+1,:) = [1 1 0]; %'y'
             end
             if ~isempty(conpsAllEchoesLvl)
                Y{end+1}   = conpsAllEchoesLvl;
                txt{end+1} = 'ConspsEchoes';
                fc(end+1,:) = [1 0 1]; % 'm'
             end   
             if ~isempty(conspsConspsEchosLvl)
                Y{end+1}   = conspsConspsEchosLvl;
                txt{end+1} = 'ConspsByConsps';
                fc(end+1,:) = [0.5 0.5 0.5]; %'grey'
             end
        case 'Conspecific'
%             histogram(ax22, conpsAllEchoesLvl, 100, 'FaceColor', 'b', 'Normalization','probability', 'DisplayName', 'All');
%             histogram(ax11, conpsAllEchoesLvl, 100, 'FaceColor', 'b', 'Normalization','probability', 'DisplayName', 'Consps Echoes');
%             histogram(ax11, conpsCallsAllLvl, 100, 'FaceColor', 'r', 'Normalization','probability', 'DisplayName', 'ConspsCalls');
%             histogram(ax11, obsAllEchoesLvl, 100, 'FaceColor', 'y', 'Normalization','probability', 'DisplayName', 'Clutter');
            Y{end+1}  = conpsAllEchoesLvl;
            txt{end+1} = 'AllConsps';
            fc(end+1,:) = [0 0 1]; % 'b'
            if ~isempty(conpsCallsAllLvl)
                Y{end+1}   = conpsCallsAllLvl;
                txt{end+1} = 'ConspsCalls';
                fc(end+1,:) = [1 0 0]; %'r';
            end
             if ~isempty(obsAllEchoesLvl)
                Y{end+1}   = obsAllEchoesLvl;
                txt{end+1} = 'Clutter';
                fc(end+1,:) = [1 1 0]; %'y'
             end    

    end % switch obsDet
% end % if
violin(Y, 'facecolor', fc, 'edgecolor', 'b', 'mc', 'k', 'medc', 'r--');
ax11.FontSize = 9;
ax11.XTick = 1:size(Y,2);
ax11.XTickLabel = txt;

title(ax11, 'Masking RxLvl')
% xlabel('dB')
ylabel('Rx lvl (dB)')

L = legend(ax11);
L.FontSize = 7;

% THE WANTED SIG
% histogram(ax11, detAllRxLvl, 100, 'FaceColor', 'g', 'Normalization','probability', 'DisplayName', 'Detected');
% histogram(ax22, detAllRxLvl, 100, 'FaceColor', 'g', 'Normalization','probability', 'DisplayName', 'Detected');

% linkaxes([ax11, ax22], 'x')

% title(ax22, 'Wanted RxLvl')
% L = legend(ax22);
% L.FontSize = 7;

% Detection Prob
ax2 = subplot(3,2,3); hold on

yyaxis left
bar(detDist, totDetAll)
ylabel('Total targets in beam')
grid minor

yyaxis right
plot(detDist, detectProb, 'LineWidth', 2);
xlabel('Detection Range [m]');
ylabel('Detection Probabilty')

title(['Detection Probabity, Beam Angle: ', num2str(rad2deg(detAngle),2),'deg'])


% polar plot
ax3 = subplot(3,2,5);
polarplot(anglesAll, distAll, '.b', 'DisplayName','all dist')
ax3 = gca;
hold on
polarplot(anglesAll(ixDetAllIn), distAll(ixDetAllIn), 'ok', 'DisplayName', 'in beam')
polarplot(detAllAngle, detAllDist, '*g', 'DisplayName', 'detections')
maxR = min([4, max(distAll)]);
rlim(ax3, [0 maxR])
% legend
title(['Polar plot, BeamRange: ', num2str(detDistIn)])

ax33 = subplot(3,2,6); hold on
scatter(detAllDist, detAllRxLvl, 18, abs(rad2deg(detAllAngle)), 'filled', 'DisplayName', 'Detected')
title('RxLvl vs Distance color by angle')
colorbar
xlabel('distance (m)')
ylabel('dB')
grid on
ax1.Position  = [0.08   0.8     0.38    0.15];
ax11.Position = [0.57  0.8     0.38    0.15];
ax2.Position  = [0.08   0.5    0.8403  0.2];
% ax22.Position = [0.5703  0.55     0.38    0.15];
ax3.Position  = [0.08  0.05   0.33   0.33];
ax33.Position = [0.57  0.06   0.33   0.35];

%%%%% Add tilte
xTitle = 1*ax1.XLim(2);
yTitle = 1.2*ax1.YLim(2);
txt = text(ax1, xTitle, yTitle, [obsDet],'FontWeight', 'bold', 'FontSize', 14);