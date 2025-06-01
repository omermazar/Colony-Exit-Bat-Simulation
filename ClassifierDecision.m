function [decisionFlg] =  ClassifierDecision(classifierResults, methods)

% decisionRule = "pks_trghs_only"; % , "All";
% decisionRule = "corrCoeffFR"; % , "All";
decisionRule = "corrCoeffMatrix"; % , "All";
switch decisionRule
    case "pks_trghs_only"
        frTH = 1.5;
        decisionFlg = classifierResults.pks_trghs_rms < frTH;
    case "corrCoeffFR"
        corrTH = 0.75;
        decisionFlg = classifierResults.corrCoeffFR >= corrTH;
    case "corrCoeffMatrix"
        corrTH = 0.8;
        decisionFlg = classifierResults.corrCoeffMatrix >= corrTH;
end % switch