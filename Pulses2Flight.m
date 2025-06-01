
function [FullVec] = Pulses2Flight(PulsesStageCell, StartPulseTimes, NumOfSamples)
%%% returns the stage in everysample of the bat from the stage in each pulse

    FullVec= repmat({''},1,NumOfSamples);

    % FullVec(StartPulseTimes) = FlightStageCell;
    for kk = 1:length(StartPulseTimes)-1
        FullVec(StartPulseTimes(kk): StartPulseTimes(kk+1)-1 ) = PulsesStageCell(kk);
    end % for kk
    FullVec(StartPulseTimes(end):end) = PulsesStageCell(end);
end % function Pulses2Flight
