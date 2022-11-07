function [loss, net] = BTMBiquadNeuralNetwork(x, trainingData, targetData, fidx, fs, nfft)

    % Create data
    if nargin < 2
        fs = 96e3;
        nfft = 4096;
        epochSize = 20e3;
        controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'saveFiles', false);
        rng(1)
        [trainingData, targetData, fvec, fNBand, fidx] = CreateBtmTrainingData(epochSize, controlparameters);
    end
    
    % Define arcitecture
    numLayers = 6;
    hiddenLayerSize = 10;
    alpha = 0.2;
    
    % Define training parameters
    miniBatchSize = 128;
    numEpochs = 100;
    
    % Loss function
    numBiquads = 2;
    lossFunc = @(net, trainingData, targetData) NNBiquadLoss(net, trainingData, targetData, numBiquads, nfft, fs, fidx, true);
    
    % Data function
    dataFunc = @() CreateBtmTrainingData(epochSize, controlparameters);

    % Train network
    [net, loss] = CreateNeuralNetwork(trainingData, targetData, numLayers, hiddenLayerSize, alpha, numEpochs, miniBatchSize, lossFunc, x);
    %[net, loss] = CreateNeuralNetwork(trainingData, targetData, numLayers, hiddenLayerSize, alpha, numEpochs, miniBatchSize, lossFunc, x, dataFunc);
end