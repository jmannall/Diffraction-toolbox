
close all
clear all

fs = 96e3;
nfft = 4096;
c = 344;

filePath = ['bayesoptResults', filesep, 'BayespotResult_Size'];
size = [700, 350, 700];
filter = {'IIR', 'Biquad', 'Biquad'};
numNetworks = length(size);

x = table('Size', [numNetworks, 5], 'VariableTypes', ["double", "double", "double", "double", "double"], 'VariableNames', ["lR", "gD", "sGD", "mG", "nL"]);
for i = 1:numNetworks
    x(i,:) = LoadBayesoptResult(size(i), filter{i});
end

% Create loss function
numFilters = 2;
nBands = 12;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);
[~, tfmag, ~, fvec, ~] = DefaultBTM(controlparameters);
[~, ~, fidx] = CreateFrequencyNBands(tfmag, fvec, nBands);

epochSize = 1e3;
numEpochs = 10;
dataFunc = @(i) CreateBtmTrainingData(epochSize, controlparameters, i);

parfor i = 1:numNetworks

    name = filter{i};
    networkSize = size(i);
    if matches(name, 'IIR')
        filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
        numOutputs = 2 * numFilters + 1;
    else
        filterFunc = @(output, target) BiquadLoss(output, target, numFilters, nfft, fs, fidx);
        numOutputs = 4 * numFilters + 1;
    end

    lossFunc = @(net, trainingData, targetData) NNFilterLoss(net, trainingData, targetData, filterFunc, true);

    [loss(i), net(i)] = CreateBTMNeuralNetwork(x(i,:), lossFunc, dataFunc, networkSize, numOutputs, numEpochs, name, false);
end

return
restarting = isfile(savePath);
if restarting
    disp('Resume training')
    load(savePath, 'net', 'losses', 'netName', 'epoch', 'index')
    nextEpoch = epoch;
    [trainingData, targetData, ~, ~, ~, index, dataPath] = CreateBtmTrainingData(epochSize, controlparameters, index);
    [numInputs, dataSize] = size(trainingData);
    numIterationsPerEpoch = floor(dataSize./miniBatchSize);
else
    disp('Start training')
    [trainingData, targetData, ~, ~, ~, index, dataPath] = CreateBtmTrainingData(epochSize, controlparameters);
    [numInputs, dataSize] = size(trainingData);
    numIterationsPerEpoch = floor(dataSize./miniBatchSize);

    nextEpoch = 1;
    losses = zeros(numEpochs, numIterationsPerEpoch, numNetworks);
    for i = 1:numNetworks
        netName(i) = {[controlparameters.filterType, '_', num2str(networkSize(i))]};
        net(i) = InitialiseNeuralNetwork(numInputs, numLayers(i), hiddenLayerSize(i), numOutputs, alpha);
    end
end

lossFunc = @(net, trainingData, targetData) NNFilterLoss(net, trainingData, targetData, filterFunc, true);

for epoch = nextEpoch:numEpochs
    disp(['Epoch ', num2str(epoch)])
    if index == 0
        [trainingData, targetData, fvec, fc, fidx, index, dataPath] = CreateBtmTrainingData(epochSize, controlparameters);
    else
        [trainingData, targetData, fvec, fc, fidx, index, dataPath] = CreateBtmTrainingData(epochSize, controlparameters, index);
    end
    save(savePath, 'net', 'losses', 'netName', 'epoch', 'index')
    for i = 1:numNetworks
    
        [net(i), losses(epoch,:, i)] = ParTrainNeuralNetwork(net(i), trainingData, targetData, epoch, miniBatchSize, numIterationsPerEpoch, lossFunc, x(i,:));
    
    end
    delete([dataPath, '.mat'])
    index = 0;
    save(savePath, 'net', 'losses', 'netName', 'epoch', 'index')
end