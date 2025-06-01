function [JamDecision] = PreyJamDecision(nTime, currPREY, AllParams)

% The function returns the Decision whether to jam or not
fs           = 1./ AllParams.SimParams.SampleTime;
nCallsTH     = AllParams.PreyFlightParams.JamBatCallsNum;
rxTH         = 10.^(AllParams.PreyFlightParams.JamSigDetectionTH/10);
timeInterval = max(AllParams.PreyFlightParams.JamIPItoReact*1e-3 * fs * (nCallsTH - 1) , AllParams.PreyFlightParams.DecisionPeriod * fs);
minIPI       = AllParams.PreyFlightParams.JamIPItoReact*1e-3 * fs;

% Jamming Criterion
relTime = max(1, nTime-timeInterval) ;

ix = find(currPREY.RxnTimes >  relTime, 1);
nTot = size(currPREY.RxnTimes,2);
if ~isempty(ix)
    ix2 = find(currPREY.RxnTimes < nTime, 1, "last");
    ixCalls = ix:ix2;
 
    rxLvls  = currPREY.RxLvls(ixCalls);
    rxTimes = currPREY.RxnTimes(ixCalls);
    detectCalls = rxLvls >= rxTH;
    
    %% the Decision Criteria
    if any(detectCalls)
        if nCallsTH == 1
            JamDecision  = true;
        else % if nCallsTH == 1
            callsIPI = diff(rxTimes(detectCalls));
            nDetections = sum(callsIPI <= minIPI);

            JamDecision = nDetections >= nCallsTH -1;
        end
    else %     if any(detectCalls)
       JamDecision  = false; 

    end %  if any(detectCalls)
    
else % ~isempty(ix)
   JamDecision  = false;
end % ~isempty(ix)




