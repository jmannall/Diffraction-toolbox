function net = InitialiseNeuralNetwork(numInputs, numLayers, hiddenLayerSize, numOutputs, alpha)

    % Create layer graph
    lgraph = layerGraph();
    
    % Create input layer
    tempLayers = [
        featureInputLayer(numInputs, 'Name', 'InputLayer')];
    lgraph = addLayers(lgraph,tempLayers);

    % Create hidden layers
    layerOut = tempLayers.Name;
    for i = 1:numLayers
        [lgraph, layerOut] = CreateHiddenLayer(lgraph, layerOut, hiddenLayerSize, alpha, i);
    end
    
    % Create output layer
    tempLayers = [
        fullyConnectedLayer(numOutputs,'Name','OutputLayer')];
    lgraph = AddConnectLayers(lgraph, tempLayers, layerOut);
    
    % clean up helper variable
    clear tempLayers;
    
    % Create network
    net = dlnetwork(lgraph);
end