function [lgraph, layerOut] = AddConnectLayers(lgraph, tempLayers, layerOut)

    lgraph = addLayers(lgraph, tempLayers);
    layerIn = tempLayers(1).Name;
    lgraph = connectLayers(lgraph,  layerOut, layerIn);
    layerOut = tempLayers(end).Name;
end