function [swarm_struct] =swarm_anlysis(BAT, AllParams, FlightInterferceSummary)


swarm_struct.num_of_detected_consp = zeros(1,AllParams.SimParams.TotalBatsNumber);
swarm_struct.num_of_masked_consp = zeros(1,AllParams.SimParams.TotalBatsNumber);
swarm_struct.num_of_AvoidBatMan = zeros(1,AllParams.SimParams.TotalBatsNumber);
swarm_struct.num_of_rel_masked_consp = zeros(1,AllParams.SimParams.TotalBatsNumber);
swarm_struct.num_of_bat_crush =  zeros(1,AllParams.SimParams.TotalBatsNumber);
swarm_struct.num_of_calls = [BAT.NumOfTimesSonarTransmits];

for kBat=1:AllParams.SimParams.TotalBatsNumber
    kCalls = swarm_struct.num_of_calls(kBat);
    swarm_struct.num_of_detected_consp(kBat) = numel([BAT(kBat).Consps_FindsStruct.DetecectedPreyWithOutInterference]);
    swarm_struct.num_of_masked_consp(kBat) = numel([BAT(kBat).Consps_FindsStruct.MaskedPreys]);
    swarm_struct.num_of_AvoidBatMan(kBat) = sum([BAT(kBat).ManueverCmdStruct(2:kCalls).ReactToBat]);
    swarm_struct.num_of_calls_rel_consps(kBat) = sum(~isnan([BAT(kBat).ManueverCmdStruct(2:kCalls).Relevant_ConspEcho_Masked]),'omitnan');
    swarm_struct.num_of_rel_masked_consp(kBat) = sum([BAT(kBat).ManueverCmdStruct(2:kCalls).Relevant_ConspEcho_Masked], 'omitnan');
    swarm_struct.num_of_bat_crush(kBat) = sum([BAT(kBat).ManueverCmdStruct(2:kCalls).BatCrush], 'omitnan');
    
end % for kBat
swarm_struct.InterferenceToDetections_bats = swarm_struct.num_of_masked_consp ./swarm_struct.num_of_detected_consp;
swarm_struct.InterferenceToRelevants_bats = swarm_struct.num_of_rel_masked_consp ./  swarm_struct.num_of_calls_rel_consps;
% swarm_struct.InterferenceToRelevants_bats = swarm_struct.num_of_rel_masked_consp./(swarm_struct.num_of_AvoidBatMan + swarm_struct.num_of_rel_masked_consp);

swarm_struct.Total_InterferenceToDetections = sum(swarm_struct.num_of_masked_consp)./sum(swarm_struct.num_of_detected_consp);
swarm_struct.Total_InterferenceToRelevants = sum(swarm_struct.num_of_rel_masked_consp)./sum(swarm_struct.num_of_calls_rel_consps);

% detection prob
detDistTest = [1,2, 3];
for ii = 1:numel(detDistTest)
    fName1 = strcat('consps', 'DetectProb_', num2str(detDistTest(ii)),'m');
    swarm_struct.(fName1) =FlightInterferceSummary.SummaryDataVectors.(fName1);
    fName2 = strcat('consps', 'TotalInBeam_', num2str(detDistTest(ii)),'m');
    swarm_struct.(fName2) =FlightInterferceSummary.SummaryDataVectors.(fName2);
end % for ii = 1:numel(detDistTest)


% swarm_struct.Total_InterferenceToRelevants = sum(swarm_struct.num_of_rel_masked_consp)./sum(swarm_struct.num_of_AvoidBatMan + swarm_struct.num_of_rel_masked_consp);
