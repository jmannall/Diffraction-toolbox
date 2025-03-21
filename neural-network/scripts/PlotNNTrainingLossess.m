function PlotNNTrainingLossess(losses, epochLosses, testLosses, figTitle)
    numIterations = length(losses);
    numEpochs = length(epochLosses);
    iterationsPerEpoch = numIterations / numEpochs;
    
    xLosses = 1:numIterations;
    xEpochLosses = iterationsPerEpoch / 2:iterationsPerEpoch:numIterations;
    figure
    plot(xLosses, losses)
    hold on
    grid on
    plot(xEpochLosses, epochLosses)
    plot(xEpochLosses, testLosses)
    xlim([0 numIterations])
    ylim([0 50])
    xlabel('Epochs')
    ylabel('Loss (dB)')
    legend('Iteration losses', 'Epoch losses', 'Test losses', 'Location', 'northeast')
    labels = split(num2str(0:20:numEpochs));
    xticks(0:20 * iterationsPerEpoch:numIterations)
    xticklabels(labels);
    title(figTitle)

    figure
    plot(xLosses, losses)
    hold on
    grid on
    plot(xEpochLosses, epochLosses)
    plot(xEpochLosses, testLosses)
    xlim([0 numIterations])
    ylim([0 10])
    xlabel('Epochs')
    ylabel('Loss (dB)')
    legend('Iteration losses', 'Epoch losses', 'Test losses', 'Location', 'northeast')
    labels = split(num2str(0:20:numEpochs));
    xticks(0:20 * iterationsPerEpoch:numIterations)
    xticklabels(labels);
    title(figTitle)
end