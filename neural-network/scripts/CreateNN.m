function [net, losses, loss] = CreateNN(hP, tP, nP)    

    net = InitialiseNN(nP.numInputs, hP.numLayers, hP.hiddenLayerSize, nP.numOutputs, nP.alpha);
    %restartFunc = @()InitialiseNeuralNetwork(numInputs, numLayers, hiddenLayerSize, numOutputs, alpha);
    
    [net, losses] = TrainNN(net, hP, tP, nP);

    loss = losses.test(end);
    disp(['Loss: ', num2str(loss)])
end