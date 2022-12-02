function UpdateNNAnimatedLinePlot(lineIterationLoss, lineEpochLoss, losses, numIterationsPerEpoch, epoch, start)

    idx = numIterationsPerEpoch * (epoch - 1) + 1:numIterationsPerEpoch * epoch;
    losses = losses(idx);
    D = duration(0,0,toc(start),'Format',"hh:mm:ss");
    epochLoss = mean(losses);
    addpoints(lineIterationLoss,idx,losses)
    addpoints(lineEpochLoss,idx(end),epochLoss)

    worker = getCurrentWorker;
    if ~isempty(worker)
        title("Worker: " + worker.ProcessId + "Epoch: " + epoch + ", Elapsed: " + string(D))
    else
        title("Epoch: " + epoch + ", Elapsed: " + string(D))
    end
    drawnow
    disp([num2str(epoch), ' Epoch loss: ', num2str(epochLoss)]);
end