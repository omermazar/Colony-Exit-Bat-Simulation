function [ ] = test_batwoman()

% This script tests the effect of the Pd on the Bats' Performance
% Nominal Variables (nBat =5, nPrey =10)
% Input- NumOfTrials : Number of trials in each set of parameters
% Output- save files of each experiment and tsummary table


% nohup matlab -nodisplay -nodesktop -nosplash -nojvm -r "cd('/home/omermazar/caveExit_Experiment'); t=CaveExit_Experiment_batwoman;exit" > log.txt &

%% Nominal Parameters load
% runPath = 'D:\Omer\BatSimulation\All Directories_backup\Bat Simulation';
% NominalParams_FilePathName = ...
%     'D:\Dropbox\University\????????\experiments\caveExit\Code\DefaultParamsTable_caveExit.xlsx';

%% for batwoman run
t = 'XXXXXXXXX Hello Experiment world XXXXXXXXXXXXX';
disp('XXXXXXXXX Hello Experiment world XXXXXXXXXXXXX')
save('helloWorld.txt', "t");