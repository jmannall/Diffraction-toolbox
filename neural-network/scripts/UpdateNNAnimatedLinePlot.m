function UpdateNNAnimatedLinePlot(lineIterationLoss, lineEpochLoss, losses, numIterationsPerEpoch, epoch, start)

    idx = numIterationsPerEpoch * (epoch - 1) + 1:numIterationsPerEpoch * epoch;
    D = duration(0, 0, toc(start), 'Format', "hh:mm:ss");
    addpoints(lineIterationLoss, idx, losses.iteration(idx))
    addpoints(lineEpochLoss, idx(end), losses.epoch(epoch))
    title("Epoch: " + epoch + ", Elapsed: " + string(D))
    drawnow

    if mod(epoch, 5) == 0
        text = ['Epoch: ' num2str(epoch), ' Epoch loss: ', num2str(losses.epoch(epoch)), ' Test loss: ', num2str(losses.test(epoch)), ' Elapsed: ', char(D)];
        worker = getCurrentWorker;
        if ~isempty(worker)
            disp(['Worker: ', num2str(worker.ProcessId), ', ', text]);
        else
            disp(text);
        end
    end
end