function [lgraph, layerOut] = CreateHiddenLayer(lgraph, layerOut, hiddenLayerSize, alpha, i)

    idx = num2str(i);
    tempLayers = [
        fullyConnectedLayer(hiddenLayerSize, 'Name', ['HiddenLayer', idx])
        batchNormalizationLayer('Name', ['NormalisationLayer_', idx])
        leakyReluLayer(alpha, 'Name', ['ActivationLayer_', idx])];
    [lgraph, layerOut] = AddConnectLayers(lgraph, tempLayers, layerOut);
end