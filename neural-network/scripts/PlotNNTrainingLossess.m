function PlotNNTrainingLossess(losses, epochLosses, figTitle)
    numIterations = length(losses);
    numEpochs = length(epochLosses);
    iterationsPerEpoch = numIterations / numEpochs;
    
    xLosses = 1:numIterations;
    xEpochLosses = 1:iterationsPerEpoch:numIterations;
    figure
    plot(xLosses, losses)
    hold on
    plot(xEpochLosses, epochLosses)
    xlim([0 numIterations])
    ylim([0 50])
    xlabel('Epochs')
    ylabel('Loss (dB)')
    legend('Iteration losses', 'Epoch losses', 'Location', 'northeast')
    labels = split(num2str(0:100:numEpochs));
    xticks(0:100 * iterationsPerEpoch:numIterations)
    xticklabels(labels);
    title(figTitle)
end