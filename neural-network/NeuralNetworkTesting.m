clc;
clear;
close all;
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultLineMarkerSize', 10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NN

% Hyperparamters
% learnRate = 1e-3;
% gradDecay = 0.9;
% sqGradDecay = 0.999;
% maxGradient = 0.9;
% 
% x.lR = learnRate;
% x.gD = gradDecay;
% x.sGD = sqGradDecay;
% x.mG = maxGradient;
% loss = XORNeuralNetwork(x);
% loss = BTMBiquadNeuralNetwork(x);

networkSize = 700; % Same as UTD
networkSize = 350;
maxLayers = networkSize / 2 ^ 2;

learnRate = optimizableVariable('lR',[1e-5 1e-1],'Type','real','Transform','log');
gradDecay = optimizableVariable('gD',[0.5 1],'Type','real');
sqGradDecay = optimizableVariable('sGD',[0.8 1],'Type','real');
maxGradient = optimizableVariable('mG',[0.5 10],'Type','real','Transform','log');
numNodes = optimizableVariable('nL',[2 min(20, maxLayers)],'Type','integer');

% Create data
fs = 96e3;
nfft = 4096;
c = 344;
epochSize = 20e3;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);

saveDir = 'runningFiles';
CheckFileDir(saveDir);
saveSeed = [saveDir, '\Seed.mat'];
saveResult = [saveDir, '\BayesoptResults.mat'];
restarting = isfile(saveSeed);
if restarting
    load(saveSeed, 'seed')
else
    seed = round(1e9 * rand(1));
    save(saveSeed, 'seed');
end
rng(seed)

[trainingData, targetData, ~, ~, ~, idx, saveData] = CreateBtmTrainingData(epochSize, controlparameters);
dataFunc = @() CreateBtmTrainingData(epochSize, controlparameters, idx);

numBiquads = 2;
nBands = 12;
[~, tfmag, ~, fvec, ~] = DefaultBTM(controlparameters);
[~, ~, fidx] = CreateFrequencyNBands(tfmag, fvec, nBands);

lossFunc = @(net, trainingData, targetData) NNBiquadLoss(net, trainingData, targetData, numBiquads, nfft, fs, fidx, true);

numOutputs = 4 * numBiquads + 1;
numEpochs = 100;
func = @(x)CreateBTMNeuralNetwork(x, lossFunc, dataFunc, networkSize, numOutputs, numEpochs);
restarting = isfile(saveResult);
disp('START')
if restarting
    load(saveResult);
    results = resume(BayesoptResults, 'SaveFileName', saveResult, 'OutputFcn', @saveToFile );
else
    result = bayesopt(func, [learnRate, gradDecay, sqGradDecay, maxGradient, numLayers], ...
    'UseParallel', true, 'MaxTime', 86400, 'MaxObjectiveEvaluations', 10, ...
    'SaveFileName', saveResult, 'OutputFcn', @saveToFile);
end

xObs = result.XAtMinObjective;
xEst = result.XAtMinEstimatedObjective;

%% Test

[lossObs, netObs] = func(xObs);
[lossEst, netEst] = func(xEst);

%% Save

save(['BayesoptResult_NNSize_', num2str(networkSize), '.mat'], "result", "xObs", "xEst", "lossObs", "netObs", "lossEst", "netEst")

%% Delete data save
delete([saveData, '.mat'])
delete(saveSeed)