function [tt] = myCateroricalTableToNumerical(tt)

oCatVars = {'ExperimentName', 'ReceiverType', 'RoomType', 'TestMode'};

catIx      = varfun(@iscategorical,tt,'output','uniform');
allCatVars = tt.Properties.VariableNames(catIx);
catVars    = setdiff(allCatVars, oCatVars);
for currCat = catVars
    tt.(currCat{1}) = double(string(tt.(currCat{1})));
end
