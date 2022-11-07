function [net, loss] = CreateNeuralNetwork(trainingData, targetData, numLayers, hiddenLayerSize, alpha, numEpochs, miniBatchSize, lossFunc, x, dataFunc)    

    [numInputs, dataSize] = size(trainingData);
    numOutputs = size(targetData, 1);

    net = InitialiseNeuralNetwork(numInputs, numLayers, hiddenLayerSize, numOutputs, alpha);

    numIterationsPerEpoch = floor(dataSize./miniBatchSize);
    
    if nargin > 9
        [net, losses] = TrainNeuralNetwork(net, trainingData, targetData, numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x, dataFunc);
    else
        [net, losses] = TrainNeuralNetwork(net, trainingData, targetData, numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x);
    end
    loss = losses(end);
end