function UpdateNNAnimatedLinePlot(lineIterationLoss, lineEpochLoss, losses, numIterationsPerEpoch, epoch, start)

    idx = numIterationsPerEpoch * (epoch - 1) + 1:numIterationsPerEpoch * epoch;
    losses = losses(idx);
    D = duration(0,0,toc(start),'Format',"hh:mm:ss");
    epochLoss = mean(losses);
    addpoints(lineIterationLoss,idx,losses)
    addpoints(lineEpochLoss,idx(end),epochLoss)

    title("Epoch: " + epoch + ", Elapsed: " + string(D))
    drawnow
    disp([num2str(epoch), ' Epoch loss: ', num2str(epochLoss)]);
end