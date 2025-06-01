%%
figure
hold on

% bat pos
% nTime = BAT.TransmittedPulsesStruct(AnalyzedPulseNum).StartPulseTime;
% xBat  = BAT.xBati(nTime);
% yBat  = BAT.yBati(nTime);
% tetaB = BAT.Teta(nTime);

% input
plot(xIn,yIn, '.k', 'MarkerSize',12)

%% clusters points
colors = lines(length(clusters));
for i = 1:length(clusters)
    plot(clusters{i}(:,1), clusters{i}(:,2), 'o', 'Color', colors(i,:), 'MarkerSize', 9) ;%, 'MarkerFaceColor', colors(i,:));
    plot(centers(i,1), centers(i,2), '*', 'Color', colors(i,:), 'MarkerSize', 9 )
end

% THe Repeatitions
plot(repeatedX, repeatedY, 'ks', 'MarkerSize', 16)
