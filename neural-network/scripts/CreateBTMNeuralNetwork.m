%% Create default neural network for BTM training

function [loss, net] = CreateBTMNeuralNetwork(x, lossFunc, dataFunc, networkSize, numOutputs, numEpochs)
    % Define arcitecture
    alpha = 0.2;
    
    % Define training parameters
    miniBatchSize = 128;

    [trainingData, targetData] = dataFunc();

    numLayers = x.nL;
    gx = 5;
    hiddenLayerSize = round((-gx + sqrt(gx ^ 2 - 4 * (-networkSize / numLayers))) / 2); % (-b + sqrt(b^2 - 4ac)) / 2a

    % Train network
    % [net, loss] = CreateNeuralNetwork(trainingData, targetData, numLayers, hiddenLayerSize, numOutputs, alpha, numEpochs, miniBatchSize, lossFunc, x, dataFunc);
    [net, loss] = CreateNeuralNetwork(trainingData, targetData, numLayers, hiddenLayerSize, numOutputs, alpha, numEpochs, miniBatchSize, lossFunc, x);
    
    idx = [num2str(x.lR), num2str(x.mG), num2str(x.gD), num2str(x.sGD), num2str(x.nL)];
    idx = erase(idx, '.');
    save(['Test_', idx, '.mat'], "net", "loss");
end