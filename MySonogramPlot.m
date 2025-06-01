
function [] = MySonogramPlot(BatDATA, DataToPlot, varargin)
% MySonogramPlot(NtimeVec, FreqVec, Powervec) will polt a singogram in time
% Imputs:
% BatDATA - a struct with the relevanr
% DataToPlot - a struct with the flags of which bat to plot
% FreqVec= the vector of the freqs in each sample
% PowerVec = the vector of the Powers
PlotSonoFlag = 1;

NumberOfBats = BatDATA.AllParams.SimParams.TotalBatsNumber;
NumberOfPreys= BatDATA.AllParams.SimParams.TotalPreysNumber;

SampleTime =  BatDATA.AllParams.SimParams.SampleTime;
MaxTime =  BatDATA.AllParams.SimParams.SimulationTime;
% MaxTime = (length(FreqVec)-1)*SampleTime;
TimeVec = [0:SampleTime:MaxTime];
DetectionTH = BatDATA.AllParams.BatSonarParams.PulseDetectionTH +...
    BatDATA.AllParams.BatSonarParams.NoiseLeveldB;


%%% PLOT OF all Bats Together
figure
clear ax1;
clear ax2;
FreqVec = zeros(NumberOfBats,length(TimeVec));
DetectionTHVec = DetectionTH*ones(1,length(TimeVec));
ax1 = cell(1, NumberOfBats);
ax2 = cell(1, NumberOfBats);%handle of the subplot axis 
for kBat = 1:NumberOfBats
       FreqVec(kBat,:) = BatDATA.BAT(kBat).BatSonarEchosMat(2,:);
       MaxFreq = max(FreqVec(kBat,:));
       MinFreq = min(nonzeros(round(FreqVec(kBat,:))));  % the min that isnt zero
       kPlot = 1+2*(kBat-1);
       ax1{1,kBat} = subplot(NumberOfBats, 2, kPlot);
       hold on
       
       % Interference from Echos
%        plot(TimeVec , 10*log10(abs(BatDATA.BAT(kBat).AllInterFromEchosPulses)), 'r.','linewidth',2)
       
       % All Interference
       plot(TimeVec , 10*log10(abs(BatDATA.BAT(kBat).AllInterPulses)), 'r-', 'linewidth',2 )
       
       % the transmitted signal and echos
       plot(TimeVec , 10*log10(abs(BatDATA.BAT(kBat).BatSonarEchosMat(1,:))),'k', 'linewidth', 2 )
       
       % Echos from the Prey
       plot(TimeVec , 10*log10(abs(BatDATA.BAT(kBat).PreyEchosVec(:))), '+m', 'linewidth',1 )
       
       % Echos from Obsticle
% % %        plot(TimeVec , 10*log10(BatDATA.BAT(kBat).ObsEchosVec(:)), '.m' )
       
       % TheDetection Thershold
       plot(TimeVec, DetectionTHVec, 'g:','linewidth',1);
       
       % Catching times of the prey
       TimesOfcatch = BatDATA.BAT(kBat).CatchPreyTimes;
       nTimesOfcatch = round(TimesOfcatch ./ SampleTime);
       plot(TimesOfcatch, 10*log10(abs(BatDATA.BAT(kBat).PreyEchosVec(nTimesOfcatch))), 'g+','MarkerSize',12,'linewidth',2)
       
       %%%  the Jammning times
       
       RefPlot= -30*ones(1,length(TimeVec));
       
%        nTimesOfInter = BatDATA.BAT(kBat).InterReportStrct.TotalInterferenceTimes;
       TimesOfInter = BatDATA.BAT(kBat).InterReportStrctOnLine.TotalInterferenceTimes;
       nTimesOfInter = round(double(TimesOfInter) ./ SampleTime);
       plot(TimesOfInter, 10*log10(abs(BatDATA.BAT(kBat).PreyEchosVec(nTimesOfInter))),'db','MarkerSize',10,'color','r','linewidth',4)
       
       title(['Bat #',num2str(kBat),' All Sonar Power'])
       axis([0 max(TimeVec) -10 120])
       
       % The Manuever Type
%        MyTimeColorByManueverPlot( RefPlot, BatDATA.BAT(kBat) ,BatDATA.AllParams, ax1{1,kBat},...
%             'Marker', 'o', 'LineStyle','none','MarkerSize',4 );
       MyTimeColorByAnyFieldPlot( RefPlot, 'FlightStageCell', BatDATA.BAT(kBat), BatDATA.AllParams, ax1{1,kBat},...
            'Marker', 'o', 'LineStyle','none','MarkerSize',4 );
       plot(TimesOfInter, RefPlot(nTimesOfInter),'db','MarkerSize',6,'color','r','linewidth',2)
       axis([0 max(TimeVec) -10 120])
       
       %%% Spectogarm
       ax2{1,kBat} = subplot(NumberOfBats, 2, kPlot+ 1);
       imagesc(10*log10(abs(BatDATA.BAT(kBat).BatSonogramWithInterference)), 'CDataMapping', 'scaled', 'xData', TimeVec);
       set(gca,'YDir','normal')
       colorbar('Location','eastoutside');
       axis([0 inf MinFreq-0.5 MaxFreq+0.5])
       title(['Bat #',num2str(kBat),' Sonogram Image'])
       Ax1(kBat) = ax1{kBat};
       Ax2(kBat) = ax2{kBat};
       
end % for kBat = 1:NumberOfBatsToPlot
linkaxes([Ax1,Ax2],'x')

