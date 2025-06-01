function [stats, callsTable] = analyzeBatalefCallTable(callsTable)

% callsTable_CaveExit_08_13_16_31
params = {'CallStartFundFreq', 'CallEndFundFreq', 'CallBW', 'Duration', 'CallPeakPower'};
callsTable.CallBW = callsTable.CallStartFundFreq - callsTable.CallEndFundFreq;

%% preliminary 
fig1 = figure; 
fig1.Position = [100 150 800 600];      
N = numel(params);
ax = cell(1,N);
txt = cell(1,N);
ncols = 2;
nrows = ceil(N/ncols);
for k=1:N
    ax{k} = subplot(nrows, ncols, k);
    stats.(params{k}).rawMean = mean(callsTable.(params{k}), 'omitnan');
    stats.(params{k}).rawStd  = std(callsTable.(params{k}), 'omitnan');
    stats.(params{k}).rawHistogram = histogram(callsTable.(params{k}),50);
    title(params{k},'Interpreter','none')
    grid on
    x1 = ax{k}.XLim(1) + 0.1*diff(ax{k}.XLim);
    y1 = ax{k}.YLim(1) + 0.9*diff(ax{k}.YLim);
    strk = [num2str(stats.(params{k}).rawMean,3), '+-' , num2str(stats.(params{k}).rawStd,3)];
    txt{k} = text(x1, y1, strk);
end
sgtitle('All Call raw')

%% Filter Outlier
ixD = callsTable.Duration > 0.002 & callsTable.Duration < 0.013;
ixFrE = callsTable.CallEndFundFreq > 20e3;
ixFrS = callsTable.CallStartFundFreq > 20e3;
ixDiff = callsTable.CallBW > 0;

ix = ixFrS & ixFrE & ixD & ixDiff;

callsTable = callsTable(ix,:);

fig2 = figure; 
fig2.Position = [500 150 800 600];

for k=1:N
    ax{k} = subplot(nrows, ncols, k);
    stats.(params{k}).mean = mean(callsTable.(params{k}), 'omitnan');
    stats.(params{k}).std  = std(callsTable.(params{k}), 'omitnan');
    stats.(params{k}).histogram = histogram(callsTable.(params{k}),50, 'FaceColor', 'g');
    title(params{k},'Interpreter','none')
    grid on
    x1 = ax{k}.XLim(1) + 0.1*diff(ax{k}.XLim);
    y1 = ax{k}.YLim(1) + 0.9*diff(ax{k}.YLim);
    strk = [num2str(stats.(params{k}).mean,3), '+-' , num2str(stats.(params{k}).std,3)];
    txt{k} = text(x1, y1, strk);
end
sgtitle('Selected Calls')