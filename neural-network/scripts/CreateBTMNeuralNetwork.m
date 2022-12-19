%% Create default neural network for BTM training

function [loss, net] = CreateBTMNeuralNetwork(x, lossFunc, dataFunc, numOutputs, numEpochs, name)
    % Define arcitecture
    alpha = 0.2;
    
    % Define training parameters
    miniBatchSize = 128;

    numLayers = x.nL;
    hiddenLayerSize = x.hL;
    gx = 5;
    numInputs = 8;

    %hiddenLayerSize = round((-gx + sqrt(gx ^ 2 - 4 * (-networkSize / numLayers))) / 2); % (-b + sqrt(b^2 - 4ac)) / 2a

    networkSize = CalculateNNIIRCost(numLayers, hiddenLayerSize, numInputs, numOutputs, gx);

    % Train network
    [trainingData, targetData] = dataFunc(1);
    [net, loss, epochLosses, losses] = CreateNeuralNetwork(trainingData, targetData, numLayers, hiddenLayerSize, numOutputs, alpha, numEpochs, miniBatchSize, lossFunc, x, dataFunc);

    idx = [num2str(x.lR), '-', num2str(x.mG), '-', num2str(x.gD), '-', num2str(x.sGD), '-', num2str(x.nL), '-', num2str(x.hL)];
    % idx = [num2str(x.lR), '-', num2str(x.nL), '-', num2str(x.nS)];
    idx = erase(idx, '.');
    CheckFileDir('NNSaves')
    save(['NNSaves', filesep, name, '-', num2str(networkSize), '_', idx, '.mat'], "net", "loss", "epochLosses", "losses");
%     save(['NNSaves', filesep, num2str(networkSize), '_', idx, '.mat'], "net", "loss", "epochLosses", "losses");
end