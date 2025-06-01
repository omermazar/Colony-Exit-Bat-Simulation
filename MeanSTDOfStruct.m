function [Stats] = MeanSTDOfStruct(StructDATA)

FNames= fieldnames(StructDATA);
for kk= 1:length(FNames)
    Name = FNames{kk};
    Field= getfield(StructDATA,FNames{kk});
    
    STRave = ['Mean',Name];
    Stats.(STRave) = mean(Field);
    
    STRstd = ['STD',Name];
    Stats.(STRstd) = std(Field);
end % for kk