%% The Refernce Table
% D:\OMER\caveExit\Figures\forFinal\PK_nBats\TableCaveExitPKnBats;
TablePk = TableCaveExitPKnBats(TableCaveExitPKnBats.MaskingByConsps == 1,:);

%% Masking Only
TablePk = TableCaveExitPKnBats(TableCaveExitPKnBats.MaskingByConsps == 1,:);
% categorical to numeric
TablePk = myCateroricalTableToNumerical(TablePk);
% Only 100 Bats 
TablePk100 = TablePk(TablePk.NumberOfBats == 100,:);

%% The original Table (1-40 Bats)
% D:\OMER\caveExit\Figures\forFinal\PK_twoEarsFA
% D:\OMER\caveExit\Figures\forFinal\PK_callLvl
% 
tt = TableCaveExit1noMask;

%% THe 100 Bats table
% D:\OMER\caveExit\Results\FalseAlarms100Bats
% D:\OMER\caveExit\Results\callLvl100BATS
% D:\OMER\caveExit\Results\velocity100Bats
tt2 = TableCaveExitPKLCaveSpeed100; % TableCaveExitPKLCaveCallLvl100; % TableCaveExit100BatsPKTwoEarswithFA;

%% find missing variables
miss12 = setdiff(tt.Properties.VariableNames, tt2.Properties.VariableNames)
miss21 = setdiff(tt2.Properties.VariableNames, tt.Properties.VariableNames)

tt2.(miss12{1}) = repmat(100, size(tt2,1),1);
tt.(miss21{1}) = repmat(0, size(tt,1),1);

%% Convert categorical to numeric
tt = myCateroricalTableToNumerical(tt);
tt2 = myCateroricalTableToNumerical(tt2);

%% Merge data-tables
tt = [tt;tt2];

%% Find missing varaibles with reference
miss12 = setdiff(tt.Properties.VariableNames, TablePk100.Properties.VariableNames)
miss21 = setdiff(TablePk100.Properties.VariableNames, tt.Properties.VariableNames)

TablePk100.TwoEarsReception =  repmat(0, size(TablePk100,1),1);
% tt2.(miss12{1}) = repmat(100, size(tt2,1),1);
% tt.(miss21{1}) = repmat(100, size(tt,1),1);

%% Merge all
tt = [tt; TablePk100];


%% Analyze the Mergerged Table
groupBy   = {'NumberOfBats', 'VelocityFree'}; % 'PulsePower'};  % 'TwoEarsReception'};
lineBy    = 'NumberOfBats'; %'PulsePower'; % 'TwoEarsReception';
subplotBy = 'none';
% preliminary
[stats] = analyzeCaveExitExp(tt, groupBy, lineBy, subplotBy, 0)

%% Final
savePath = 'D:\OMER\caveExit\FinalFIgures\Speed';
[stats] = analyzeCaveExitExp(tt, groupBy, lineBy, subplotBy, 1, savePath)