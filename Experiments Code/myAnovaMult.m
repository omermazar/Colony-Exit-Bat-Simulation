% myMultiAnova_CallLvl
tt = TableCaveExitcallLvl(TableCaveExitcallLvl.NumberOfBats == 100,:);
[p,tbl,sAnova] = anova1(tt.CrushesObsTotal, tt.PulsePower);
[results, m, h, gNames] = multcompare(sAnova);
tblResults = array2table(results, "VariableNames", ["A", "B", "LowerLimit", "AminusB", "UpperLimit", "pVal"]);
