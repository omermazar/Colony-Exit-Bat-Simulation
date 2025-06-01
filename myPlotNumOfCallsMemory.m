SampleTime = BatDATA.AllParams.SimParams.SampleTime;
mem = 0.75;

figure

subplot(2,1,1); hold on; grid on
nCalls = cell(size(BatDATA.BAT, 2),1);

for batNum = 1:size(BatDATA.BAT, 2)

    numCalls = BatDATA.BAT(batNum).NumOfTimesSonarTransmits;
    memSamples = mem/SampleTime;
    allTimes = vertcat(BatDATA.BAT(batNum).TransmittedPulsesStruct.StartPulseTime);
    nCalls{batNum}   = zeros(size(allTimes));
    for k = 1:numCalls
        currT = allTimes(k);
        kStart = max(0, currT-memSamples);
        nCalls{batNum}(k) = sum(allTimes <= currT & allTimes >= kStart);
    end

    plot(allTimes*SampleTime, nCalls{batNum},'.-')

end
ylabel('calls in mem')
xlabel('time (sec)')
title('calls im Memory per Bat')


subplot(2,1,2)
histogram(vertcat(nCalls{:}),100)
ylabel('accurances')
xlabel('calls in mem')
title('histogram of all Bats')
sgtitle(strcat('Memory = ', num2str(mem), 'seconds'))
