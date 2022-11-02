clc;
clear;
close all;
set(0, 'DefaultLineLineWidth', 2);
set(groot, 'defaultLineMarkerSize', 10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NN

% Hyperparamters
% learnRate = 1e-2;
% gradDecay = 0.9;
% sqGradDecay = 0.999;
% maxGradient = 5;

% x.lR = learnRate;
% x.gD = gradDecay;
% x.sGD = sqGradDecay;
% x.mG = maxGradient;
% loss = XORNeuralNetwork(x);

learnRate = optimizableVariable('lR',[1e-5 1e-1],'Type','real','Transform','log');
gradDecay = optimizableVariable('gD',[0.5 1],'Type','real');
sqGradDecay = optimizableVariable('sGD',[0.8 1],'Type','real');
maxGradient = optimizableVariable('mG',[0.5 10],'Type','real','Transform','log');

func = @(x)XORNeuralNetwork(x);
result = bayesopt(func, [learnRate, gradDecay, sqGradDecay, maxGradient], 'UseParallel', true);

xObs = result.XAtMinObjective;
xEst = result.XAtMinEstimatedObjective;

[lossObs, netObs] = XORNeuralNetwork(xObs);
[lossEst, netEst] = XORNeuralNetwork(xEst);
