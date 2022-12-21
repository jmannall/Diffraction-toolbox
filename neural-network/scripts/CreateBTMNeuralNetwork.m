%% Create default neural network for BTM training

function [loss, net] = CreateBTMNeuralNetwork(x, lossFunc, dataFunc, numOutputs, numEpochs, saveFile)
    % Define arcitecture
    alpha = 0.2;
    
    % Define training parameters
    miniBatchSize = 128;

    % Train network
    [trainingData, targetData] = dataFunc(1);
    [net, loss, epochLosses, losses] = CreateNeuralNetwork(trainingData, targetData, x.nL, x.hL, numOutputs, alpha, numEpochs, miniBatchSize, lossFunc, x, dataFunc);

    CheckFileDir('NNSaves')
    save(['NNSaves', filesep, saveFile, '.mat'], "net", "loss", "epochLosses", "losses");
end