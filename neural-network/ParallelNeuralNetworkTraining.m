
close all
clear all

fs = 96e3;
nfft = 16384;
c = 344;
filePath = ['bayesoptResults', filesep, 'BayespotResult_Size'];
% size = [350, 700, 350, 700];
% filter = {'IIR', 'IIR', 'Biquad', 'Biquad'};

size = [200; 1e3; 5e3; 10e3; 20e3;];
size = [size; size; size; size];

iir = 'IIR';
biquad = 'Biquad';
filter = [{iir; iir; iir; iir; iir; biquad; biquad; biquad; biquad; biquad}];
filter = [filter; filter];
numNetworks = length(size);

% x = table('Size', [numNetworks, 5], 'VariableTypes', ["double", "double", "double", "double", "double"], 'VariableNames', ["lR", "gD", "sGD", "mG", "nL"]);
% for i = 1:numNetworks
%     x(i,:) = LoadBayesoptResult(size(i), filter{i});
% end


lR(1:20,1) = 1e-3;
gD(1:20,1) = 0.9;
sGD(1:20,1) = 0.99;
mG(1:20,1) = 1;
nL(1:10,1) = 2;
nL(11:20,1) = 5;

x = table(lR, gD, sGD, mG, nL, 'VariableNames', ["lR", "gD", "sGD", "mG", "nL"]);

% Create loss function
numFilters = 2;
nBands = 12;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);
[~, tfmag, ~, fvec, ~] = DefaultBTM(controlparameters);
[~, ~, fidx] = CreateFrequencyNBands(tfmag, fvec, nBands);

epochSize = 20e3;
numEpochs = 500;
dataFunc = @(i) CreateBtmTrainingData(epochSize, controlparameters, i);

disp('Start parallel training')
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

% restarting = isfile(savePath);
% if restarting
%     disp('Resume training')
%     load(savePath, 'net', 'losses', 'netName', 'epoch', 'index')
%     nextEpoch = epoch;
%     [trainingData, targetData, ~, ~, ~, index, dataPath] = CreateBtmTrainingData(epochSize, controlparameters, index);
%     [numInputs, dataSize] = size(trainingData);
%     numIterationsPerEpoch = floor(dataSize./miniBatchSize);
% else
%     disp('Start training')
%     [trainingData, targetData, ~, ~, ~, index, dataPath] = CreateBtmTrainingData(epochSize, controlparameters);
%     [numInputs, dataSize] = size(trainingData);
%     numIterationsPerEpoch = floor(dataSize./miniBatchSize);
% 
%     nextEpoch = 1;
%     losses = zeros(numEpochs, numIterationsPerEpoch, numNetworks);
%     for i = 1:numNetworks
%         netName(i) = {[controlparameters.filterType, '_', num2str(networkSize(i))]};
%         net(i) = InitialiseNeuralNetwork(numInputs, numLayers(i), hiddenLayerSize(i), numOutputs, alpha);
%     end
% end
% 
% lossFunc = @(net, trainingData, targetData) NNFilterLoss(net, trainingData, targetData, filterFunc, true);
% 
% for epoch = nextEpoch:numEpochs
%     disp(['Epoch ', num2str(epoch)])
%     if index == 0
%         [trainingData, targetData, fvec, fc, fidx, index, dataPath] = CreateBtmTrainingData(epochSize, controlparameters);
%     else
%         [trainingData, targetData, fvec, fc, fidx, index, dataPath] = CreateBtmTrainingData(epochSize, controlparameters, index);
%     end
%     save(savePath, 'net', 'losses', 'netName', 'epoch', 'index')
%     for i = 1:numNetworks
%     
%         [net(i), losses(epoch,:, i)] = ParTrainNeuralNetwork(net(i), trainingData, targetData, epoch, miniBatchSize, numIterationsPerEpoch, lossFunc, x(i,:));
%     
%     end
%     delete([dataPath, '.mat'])
%     index = 0;
%     save(savePath, 'net', 'losses', 'netName', 'epoch', 'index')
% end