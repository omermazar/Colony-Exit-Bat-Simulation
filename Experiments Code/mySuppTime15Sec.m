function [stats] = mySuppTime15Sec(TableCaveExit)

% input:
% C:\Users\YossiYNB3\Documents\University\Research\experiments\caveExit\Figures\Final01\PK_RM_Comp\TableCaveExitPkRm.txt
groupBy = {'NumberOfBats', 'BatSpecies' };
lineBy = 'NumberOfBats';
checkParams = {'obsTotalInterferenceRatio', 'dist2Consps_Mean'};
uLines = [5, 10, 40, 100];
cLines = lines(numel(uLines));

%% filter in only PK
TableCaveExit = TableCaveExit(TableCaveExit.BatSpecies == 'PK',:);

%% Grouping
[g, TID]= findgroups(TableCaveExit(:, groupBy));
TableCaveExit.group = g;
analysisTable = TID;
analysisTable.numel  = splitapply(@numel, TableCaveExit.NumberOfBats, g);

G = unique(g);
Y = cell(1,numel(G));
xl = cell(1,numel(G));
cLines = lines(numel(G));

for kGroup = G'
    xl{kGroup} = num2str(unique(TableCaveExit.NumberOfBats(g == G(kGroup))));
    Y{kGroup} = TableCaveExit.ExitTimesSec(g == G(kGroup));
end % for
figure;

violin(Y, 'facecolor',cLines, 'edgecolor','b',...
'bw',0.3,...
'mc','k',...
'medc','r:')
grid on
xticklabels(xl)
ylabel('exit-time (sec)')
xlabel('Number of Bats')
legend off
