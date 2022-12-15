clear all
close all

load(['NNSaves', filesep, 'IIR-10000_0001-1-09-099-5.mat'])

fs = 96e3;
nfft = 16384;
c = 344;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);

epochSize = 20e3;
[trainingData, targetData, ~, ~, fidx] = CreateBtmTrainingData(epochSize, controlparameters, epochSize);

dataFunc = @(idx)CreateBtmTrainingData(epochSize, controlparameters, idx);
numIIRFilters = 2;
gamma = 4;
filterFunc = @(output, target)IIRFilterLoss(output, target, numIIRFilters, nfft, fs, fidx);
NNLossFunc = @(net, trainingData, targetData)NNFilterLoss(net, trainingData, targetData, filterFunc, true);
lossFunc = @(net, psi, psiInv, trainingData, targetData)TransformationLoss(net, psi, psiInv, trainingData, targetData, NNLossFunc, numIIRFilters, gamma);
miniBatchSize = 128;
numEpochs = 10;
numIterationsPerEpoch = floor(epochSize./miniBatchSize);

[psi, epochLosses, losses] = TrainTransformation(net, trainingData, targetData, numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, dataFunc);
