
close all
clear all

%% Data

fs = 48e3;
nfft = 8192;
c = 344;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);

filePath = ['bayesoptResults', filesep, 'BayespotResult_Size'];
saveDir = 'NNSaves';
CheckFileDir(saveDir)

iir = '2_iir';
iirW = '2_iirW';

%% NN complexity

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

%% Create input variables

numSizes = length(networkSizes);
numLayers = length(layers);
numNetworks = numLayers * numSizes;

lR(1:numNetworks,1) = 1e-3;
gD(1:numNetworks,1) = 0.9;
sGD(1:numNetworks,1) = 0.99;
mG(1:numNetworks,1) = 1;
nL = layers .* ones(numSizes, numLayers);
nL = reshape(nL, numNetworks, 1);
hL = reshape(hL, numNetworks, 1);

x = table(lR, gD, sGD, mG, nL, hL, 'VariableNames', ["lR", "gD", "sGD", "mG", "nL", "hL"]);
x = [x; x];

%% Create loss function

numFilters = 2;
nBands = 8;
[~, tfmag, ~, fvec, ~] = DefaultBTM(controlparameters);
[~, fc, fidx] = CreateFrequencyNBands(tfmag, fvec, nBands);
controlparameters.fs = 2 * fs;
controlparameters.nfft = 2 * nfft;

epochSize = 20e3;
numEpochs = 500;

filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
numOutputs = 2 * numFilters + 1;
weight = 20;

%% Run grid search

rng(2)
disp('Start parallel training')
parfor i = 1:2 * numNetworks
    
    if i > numNetworks
        dataFunc = @(i) CreateBtmTrainingDataWeighted(epochSize, controlparameters, weight, i);
        name = iirW;
    else        
        dataFunc = @(i) CreateBtmTrainingData(epochSize, controlparameters, i);
        name = iir;
    end
    input = x(i,:);

    gx = 5;
    numInputs = 8;
    networkSize = CalculateNNIIRCost(input.nL, input.hL, numInputs, numOutputs, gx);
    %hiddenLayerSize = round((-gx + sqrt(gx ^ 2 - 4 * (-networkSize / numLayers))) / 2); % (-b + sqrt(b^2 - 4ac)) / 2a

    idx = [num2str(input.lR), '-', num2str(input.mG), '-', num2str(input.gD), '-', num2str(input.sGD), '-', num2str(input.nL), '-', num2str(input.hL)];
    idx = erase(idx, '.');
    saveFile = [name, '-', num2str(networkSize), '_', idx];
    networkExists = exist([cd filesep saveDir filesep saveFile '.mat'], "file");
    if networkExists == 2
        disp('Already trained')
    else
        worker = getCurrentWorker;
        disp(['Begin worker: ' num2str(worker.ProcessId)])
    
        lossFunc = @(net, trainingData, targetData) NNFilterLoss(net, trainingData, targetData, filterFunc, true);
    
        [lossStore(i), netStore(i)] = CreateBTMNeuralNetwork(input, lossFunc, dataFunc, numOutputs, numEpochs, saveFile);
    end
end

saveDir = 'gridSearchResults';
CheckFileDir(saveDir)

loss = lossStore(1:numNetworks);
net = netStore(1:numNetworks);
save([saveDir, filesep, 'NN-', num2str(round(fs / 1e3))], "loss", "net")

loss = lossStore(numNetworks + 1:end);
net = netStore(numNetworks + 1:end);
save([saveDir, filesep, 'NN-', num2str(round(fs / 1e3)), '_Weighted'], "loss", "net")
