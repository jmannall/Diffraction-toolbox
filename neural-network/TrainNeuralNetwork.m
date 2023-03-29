function TrainNeuralNetwork(numLayers, size, learnRate, saveDir)
    
    % Control parameters
    fs = 48e3;
    nfft = 8192;
    c = 344;
    
    controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2, 'noDirect', true);
    
    % Filter parameters
    numFilters = 2;
    nBands = 8;
    [~, tfmag, ~, fvec, ~] = DefaultBTM(controlparameters);
    [~, ~, fidx] = CreateFrequencyNBands(tfmag, fvec, nBands);
    
    filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
    
    controlparameters.fs = 2 * fs;
    controlparameters.nfft = 2 * nfft;
    
    % Training parameters
    numEpochs = 5;
    epochSize = 20e3;
    miniBatchSize = 128;

    gradDecay = 0.9;
    sqGradDecay = 0.99;
    maxGrad = 0.9;

    lossFunc = @(net, trainingData, targetData) NNFilterLoss(net, trainingData, targetData, filterFunc, true);
    dataFunc = @(idx) CreateBtmTrainingData(epochSize, controlparameters, idx);
    testFunc = @(net, trainingData, targetData) NNFilterLoss(net, trainingData, targetData, filterFunc, false);
    
    % Network parameters
    numInputs = 8;
    numOutputs = 2 * numFilters + 1;
    alpha = 0.2;
    
    saveFile = [name, '-', num2str(networkSize), '_', idx];
    savePath = [saveDir filesep saveFile];

    % NN complexity
    gx = 5;
    CcPZTransform = 6;
    numPZ = numOutputs - 1;
    CcKTransform = 1;
    CcOutputTransform = numPZ * CcPZTransform + CcKTransform;
    CcIIR = 6;
    
    a = numLayers - 1;
    b = numLayers + numLayers * gx + numInputs + numOutputs;
    c = CcIIR + CcOutputTransform + numOutputs - size;
    
    hiddenLayerSize = round((-b + sqrt(b .^ 2 - 4 .* c * a)) ./ (2 * a));
    
    % Parameter structs
    hyperParameters = struct('learnRate', learnRate, 'numLayers', numLayers, 'hiddenLayerSize', hiddenLayerSize);
    trainingParameters = struct('lossFunc', lossFunc, 'dataFunc', dataFunc, 'testFunc', testFunc, 'numEpochs', numEpochs, 'epochSize', epochSize, ...
        'miniBatchSize', miniBatchSize, 'gradDecay', gradDecay, 'sqGradDecay', sqGradDecay, 'maxGrad', maxGrad);
    networkParamters = struct('numInputs', numInputs, 'numOutputs', numOutputs, 'alpha', alpha, 'savePath', savePath);
    
    % Train the network
    [loss, net] = CreateNN(hyperParameters, trainingParameters, networkParamters);
end
