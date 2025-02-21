function [net, losses, loss] = CreateNN(hP, tP, nP, net, losses) 

    numIterationsPerEpoch = floor(tP.epochSize./tP.miniBatchSize);
    if (nargin < 4)
        net = InitialiseNN(nP.numInputs, hP.numLayers, hP.hiddenLayerSize, nP.numOutputs, nP.alpha);
        %restartFunc = @()InitialiseNeuralNetwork(numInputs, numLayers, hiddenLayerSize, numOutputs, alpha);
        losses.iteration = zeros(1, tP.numEpochs * numIterationsPerEpoch);
        [losses.epoch, losses.test] = deal(1e3 * ones(1, tP.numEpochs));
    else
        numIterations = length(losses.iteration);
        losses.iteration = [losses.iteration, zeros(1, tP.numEpochs * numIterationsPerEpoch - numIterations)];
        numEpochs = length(losses.epoch);
        losses.epoch = [losses.epoch, deal(1e3 * ones(1, tP.numEpochs - numEpochs))];
        losses.test = [losses.test, deal(1e3 * ones(1, tP.numEpochs - numEpochs))];
    end

    [net, losses] = TrainNN(net, hP, tP, nP, losses);

    loss = losses.test(end);
    disp(['Loss: ', num2str(loss)])
end