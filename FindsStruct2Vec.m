function [xFinds,yFinds,nTimesFinds] = FindsStruct2Vec(NumOfFinds, FindsStruct, xyResolution, ObsOrPreyFlag)
% this function creates vectors ofx-y  findings from the struct of the
% findings
% Inputs:
%       NumOfFinds - the Toal times the bat found target
%       FindsStruct - the struct of the finnings with the fields-
%       xFinds, yFinds, nTimesFinds, TotalNumOfTimesFinds

%%% NumOfFinds = FindsStruct.TotalNumOfTimesFinds;
% try
%     FindsStruct.nTimesFinds = FindsStruct.DetectedTimes;
% % catch
%     
% end
switch ObsOrPreyFlag
    case 'Prey'
        Pulses = [FindsStruct(find([FindsStruct.IsAnyPreyDetected])).TransmittedPulseNum];
        xFinds = [FindsStruct(find([FindsStruct.IsAnyPreyDetected])).xFinds] * xyResolution ;
        yFinds = [FindsStruct(find([FindsStruct.IsAnyPreyDetected])).yFinds] * xyResolution;
        nTimesFinds = [FindsStruct.DetectedTimes]; 
   
    case 'Obs' % to be recoded !!!! but DO the JOB
        nn=1;
        nTimesFinds = zeros(1,NumOfFinds);
        xFinds = zeros(1,NumOfFinds*10);
        yFinds = xFinds;
        for k = 1:NumOfFinds
            NumberOfCurrentFinds = length(FindsStruct(k).xFinds);         
            if NumberOfCurrentFinds > 0
                xFinds(nn:(nn+NumberOfCurrentFinds-1)) = FindsStruct(k).xFinds .* xyResolution;
                yFinds(nn:(nn+NumberOfCurrentFinds-1)) = FindsStruct(k).yFinds .* xyResolution;
%                 nTimesFinds(nn:(nn+NumberOfCurrentFinds-1)) = FindsStruct(k).nTimesFinds;
            end % if NumberOfCurrentFinds
            nn = nn+ NumberOfCurrentFinds;
        end % for k
        xFinds = nonzeros(xFinds);
        yFinds = nonzeros(yFinds);
end % Switch