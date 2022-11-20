
close all
clear all

fs = 96e3;
nfft = 4096;
c = 344;
epochSize = 20e3;
miniBatchSize = 256;
numEpochs = 500;
alpha = 0.2;

load ('BayesoptResult_Size_700_Filter_IIR.mat', 'xEst', 'xObs');
x = [xEst; xObs];

numNetworks = size(x, 1);
networkSize = [700, 700];
numLayers = x.nL;
gx = 5;
hiddenLayerSize = round((-gx + sqrt(gx ^ 2 - 4 * (-networkSize ./ numLayers))) / 2); % (-b + sqrt(b^2 - 4ac)) / 2a

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);
CheckFileDir('data')
savePath = 'data\NeuralNetworkTest.mat';

% Create loss function
numFilters = 2;
nBands = 12;
[~, tfmag, ~, fvec, ~] = DefaultBTM(controlparameters);
[~, ~, fidx] = CreateFrequencyNBands(tfmag, fvec, nBands);
biquad = true;
if biquad
    filterFunc = @(output, target) BiquadLoss(output, target, numFilters, nfft, fs, fidx);
    numOutputs = 4 * numFilters + 1;
    controlparameters.filterType = 'Biquad';
else
    filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
    numOutputs = 2 * numFilters + 1;
    controlparameters.filterType = 'IIR';
end

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