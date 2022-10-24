function [net, loss] = CreateNeuralNetwork(trainingData, targetData, numLayers, hiddenLayerSize, alpha, numEpochs, miniBatchSize, lossFunc, x)    

    [numInputs, dataSize] = size(trainingData);
    numOutputs = size(targetData, 1);

    net = InitialiseNeuralNetwork(numInputs, numLayers, hiddenLayerSize, numOutputs, alpha);

    numIterationsPerEpoch = floor(dataSize./miniBatchSize);
    
    [net, losses] = TrainNeuralNetwork(net, trainingData, targetData, numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x);
    loss = losses(end);
end