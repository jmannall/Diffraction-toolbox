function [net, loss, epochLosses, losses] = CreateNeuralNetwork(trainingData, targetData, numLayers, hiddenLayerSize, numOutputs, alpha, numEpochs, miniBatchSize, lossFunc, x, dataFunc)    

    [numInputs, dataSize] = size(trainingData);

    net = InitialiseNeuralNetwork(numInputs, numLayers, hiddenLayerSize, numOutputs, alpha);
    restartFunc = @()InitialiseNeuralNetwork(numInputs, numLayers, hiddenLayerSize, numOutputs, alpha);

    numIterationsPerEpoch = floor(dataSize./miniBatchSize);
    
    if nargin > 10
        [net, epochLosses, losses] = TrainNeuralNetwork(net, trainingData, targetData, numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x, restartFunc, dataFunc);
    else
        [net, epochLosses, losses] = TrainNeuralNetwork(net, trainingData, targetData, numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x, restartFunc);
    end
    disp(epochLosses)
    loss = epochLosses(end);
end