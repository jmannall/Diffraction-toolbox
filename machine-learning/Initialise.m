%% Initialise custom neural network

function output = Initialise(trainingData, targetData, learnRate, numEpochs, miniBatchSize, numInput, numHidden, numOutput, Activation, Delactivation, range)
    numLayers = length(numHidden);
    numIterationsPerEpoch = floor(size(trainingData,1)./miniBatchSize);

    [weightsLayer, biasLayer] = deal(cell(1, numLayers));
    weightsLayer{1} = -range / 2 + range * rand(numInput, numHidden(1));
    biasLayer{1} = -range / 2 + range * rand(1, numHidden(1));
    for i = 2:numLayers
        weightsLayer{i} = -range / 2 + range * rand(numHidden(i - 1), numHidden(i));
        biasLayer{i} = -range / 2 + range * rand(1, numHidden(i));
    end
    weightsOutput = -range / 2 + range * rand(numHidden(numLayers), numOutput);
    biasOutput = -range / 2 + range * rand(1, numOutput);

    output = struct('trainingData', trainingData, 'targetData', targetData, 'learnRate', learnRate, 'numEpochs', numEpochs, 'miniBatchSize', miniBatchSize, 'numIterationsPerEpoch', numIterationsPerEpoch, 'numLayers', numLayers, ...
    'weightsLayer', {weightsLayer}, 'weightsOutput', weightsOutput, 'biasLayer', {biasLayer}, 'biasOutput', biasOutput, 'losses', [], ...
    'hidden', [], 'hiddenOut', [], 'Activation', Activation, 'Delactivation', Delactivation, 'output', [], 'outputFinal', []);
end