function tableAll = concatenateMaskingconfusionTables(tableMasking, tableConf)
% concatenate Masking No-Masking and Confusion -cluster Tables

%% define the fields to concatenate
reqFields = ["NumberOfBats", "MaskingByConsps", "ConfusionObs", "ObsClusterFlag", "ExpTypeConfCluster" "ExitSuccess", "ExitTimesSec", "CrushesObsPerSec"];

%% tableConf
% Check if field exisists
if ~ismember('CrushesObsPerSec', tableConf.Properties.VariableNames)
    totalFlightTime = tableConf.ExitTimesSec;
    totalFlightTime(isnan(tableConf.ExitTimesSec)) = 15;
    tableConf.CrushesObsPerSec = tableConf.CrushesObsTotal ./ totalFlightTime;
end
% from the confsion-cluster table take all data from the required fields
tableAll = tableConf(:,reqFields);

%% tableMasking
% rows without masking
ixNoMasking = tableMasking.MaskingByConsps == 0; 

% Fill Missing DATA
% Check if field exsists
if ~ismember('CrushesObsPerSec', tableMasking.Properties.VariableNames)
    totalFlightTime = tableMasking.ExitTimesSec;
    totalFlightTime(isnan(tableMasking.ExitTimesSec)) = 15;
    tableMasking.CrushesObsPerSec = tableMasking.CrushesObsTotal ./ totalFlightTime;
end
hMasking= height(tableMasking);
tableMasking.ConfusionObs       = false(hMasking,1);
tableMasking.ObsClusterFlag     = false(hMasking,1);
tableMasking.ExpTypeConfCluster = 0*ones(hMasking,1);
tableMasking.ExpTypeConfCluster(ixNoMasking) = -1*ones(sum(ixNoMasking),1);

% add all data with no masking
ixNoMasking = tableMasking.MaskingByConsps == 0; 
tempMasking = tableMasking(ixNoMasking, reqFields);

% add one-hundred Bats with Masking and no confusion to tableAll
ix100Masking = ~ixNoMasking  & tableMasking.NumberOfBats == 100;
tempMasking = [tempMasking; tableMasking(ix100Masking, reqFields)];

%% Conc both tables
tableAll = [tableAll; tempMasking];