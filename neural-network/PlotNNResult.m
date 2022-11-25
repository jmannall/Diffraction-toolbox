

numIterations = length(losses);
numEpochs = length(epochLosses);
iterationsPerEpoch = numIterations / numEpochs;

xLosses = 1:numIterations;
xEpochLosses = 1:iterationsPerEpoch:numIterations;
figure
plot(xLosses, losses)
hold on
plot(xEpochLosses, epochLosses)