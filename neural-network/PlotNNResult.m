
%close all

numIterations = length(losses);
numEpochs = length(epochLosses);
iterationsPerEpoch = numIterations / numEpochs;

xLosses = 1:numIterations;
xEpochLosses = 1:iterationsPerEpoch:numIterations;
figure
plot(xLosses, losses)
hold on
plot(xEpochLosses, epochLosses)

cost = [];
runningMean = [];
for i = 1:numEpochs
    idx = max(i - 5, 1):i - 1;
    runningMean(i) = mean(epochLosses(idx));
    cost(i) = runningMean(i) - epochLosses(i);
end

xEpoch = 1:numEpochs;
figure
plot(xEpoch, runningMean)

figure
plot(xEpoch, cost)
