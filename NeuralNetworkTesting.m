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

learnRate = optimizableVariable('lR',[1e-5 1e-1],'Type','real','Transform','log');
gradDecay = optimizableVariable('gD',[0.5 1],'Type','real');
sqGradDecay = optimizableVariable('sGD',[0.8 1],'Type','real');
maxGradient = optimizableVariable('mG',[0.5 10],'Type','real','Transform','log');

% Create data
fs = 96e3;
nfft = 4096;
epochSize = 20e3;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'saveFiles', false);
rng(1)
[trainingData, targetData, fvec, fNBand, fidx] = CreateBtmTrainingData(epochSize, controlparameters);

func = @(x)BTMBiquadNeuralNetwork(x, trainingData, targetData, fidx, fs, nfft);
result = bayesopt(func, [learnRate, gradDecay, sqGradDecay, maxGradient], 'UseParallel', true, 'MaxTime', 72000, 'MaxObjectiveEvaluations', 30);

xObs = result.XAtMinObjective;
xEst = result.XAtMinEstimatedObjective;

%% Save

save('NNTestRun', "result", "xObs", "xEst", "netObs", "netEst")

%% Test

[lossObs, netObs] = BTMBiquadNeuralNetwork(xObs);
[lossEst, netEst] = BTMBiquadNeuralNetwork(xEst);
