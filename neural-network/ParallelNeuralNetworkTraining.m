
close all
clear all

fs = 48e3;
nfft = 8192;
c = 344;
controlparameters = struct('fs', 2 * fs, 'nfft', 2 * nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);

filePath = ['bayesoptResults', filesep, 'BayespotResult_Size'];

iir = 'IIR';
iirW = 'iirW';

nI = 8;
nO = 5;
gx = 5;

networkSizes = [2e3; 4e3; 6e3; 8e3; 10e3];
layers = [2, 3, 4, 5, 6, 7];

CcPZTransform = 5;
numPZ = nO - 1;
CcKTransform = 1;
CcOutputTransform = numPZ * CcPZTransform + CcKTransform;

CcIIR = 6;

a = layers - 1;
b = layers + layers * gx + nI + nO;
c = CcIIR + CcOutputTransform + nO - networkSizes;
hL = round((-b + sqrt(b .^ 2 - 4 .* c * a)) ./ (2 * a));

% const = ones(4, 1);
% nL(1:20,1) = [2 * const; 3 * const; 5 * const; 6 * const; 7 * const];
% layers1 = [32; 48; 64; 80];
% layers2 = [20; 32; 44; 56];
% layers3 = [16; 24; 32; 40];
% layers4 = [12; 20; 28; 36];
% layers5 = [10; 18; 24; 32];
% hL(1:20,1) = [layers1; layers2; layers3; layers4; layers5];

const = ones(20,1);
numSizes = length(networkSizes);
numLayers = length(layers);
numNetworks = numLayers * numSizes;
for i = 1:numSizes
    for j = 1:numLayers
        cost(i, j) = CalculateNNIIRCost(layers(j), hL(i, j), nI, nO, gx);
    end
end
x = 1:numSizes;
figure
plot(x, cost)
legend('2', '3', '4', '5', '6', '7')

lR(1:numNetworks,1) = 1e-3;
gD(1:numNetworks,1) = 0.9;
sGD(1:numNetworks,1) = 0.99;
mG(1:numNetworks,1) = 1;
nL = layers .* ones(numSizes, numLayers);
nL = reshape(nL, numNetworks, 1);
hL = reshape(hL, numNetworks, 1);

x = table(lR, gD, sGD, mG, nL, hL, 'VariableNames', ["lR", "gD", "sGD", "mG", "nL", "hL"]);
x = [x; x];

% Create loss function
numFilters = 2;
nBands = 8;
[~, tfmag, ~, fvec, ~] = DefaultBTM(controlparameters);
[~, ~, fidx] = CreateFrequencyNBands(tfmag, fvec, nBands);

epochSize = 20e3;
numEpochs = 500;

filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
numOutputs = 2 * numFilters + 1;
weight = 20;

disp('Start parallel training')
parfor i = 1:2 * numNetworks

    if i > numNetworks
        dataFunc = @(i) CreateBtmTrainingDataWeighted(epochSize, controlparameters, weight, i);
        name = iir;
    else        
        dataFunc = @(i) CreateBtmTrainingData(epochSize, controlparameters, i);
        name = iirW;
    end

    worker = getCurrentWorker;
    disp(['Begin worker: ' num2str(worker.ProcessId)])

    lossFunc = @(net, trainingData, targetData) NNFilterLoss(net, trainingData, targetData, filterFunc, true);

    [loss(i), net(i)] = CreateBTMNeuralNetwork(x(i,:), lossFunc, dataFunc, numOutputs, numEpochs, name);
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