function [loss, net] = XORNeuralNetwork(x)

    % Create data
    trainingData = ([0 0
        1 0
        0 1
        1 1])';
    targetData = [1 0 0 1];
    
    % Define arcitecture
    numLayers = 1;
    hiddenLayerSize = 2;
    alpha = 0.2;
    
    % Define training parameters
    miniBatchSize = 4;
    numEpochs = 100;
    
    % Loss function
    lossFunc = @(net, trainingData, targetData) MeanSquaredError(net, trainingData, targetData);
    
    % Train network
    [net, loss] = CreateNeuralNetwork(trainingData, targetData, numLayers, hiddenLayerSize, alpha, numEpochs, miniBatchSize, lossFunc, x);
end