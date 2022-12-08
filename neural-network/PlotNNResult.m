
close all
% clear all

%% Load networks

loadDir = 'NNSaves';
listing = dir(loadDir);
listing = {listing(3:end).name};
networks.filterType = extractBefore(listing, '-');
networks.networkSize = extractBetween(listing, '-', '_');
networks.numLayers = extractBetween(listing, '099-', '.');

numNetworks = length(listing);

%% Create Data

fs = 96e3;
nfft = 16384;
c = 344;
testSize = 1e3;

controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);
[inputData, targetData] = CreateBtmTrainingData(testSize, controlparameters, testSize);
targetData = extractdata(targetData);

numFilters = 2;
nBands = 12;
[~, tfmag, ~, fvec, ~] = DefaultBTM(controlparameters);
[~, fc, fidx] = CreateFrequencyNBands(tfmag, fvec, nBands);
numFreq = fidx(end);

%% Plots

close all

numSelect = 5;
colours = colororder;
colours = colours(1:numSelect, :);
lossFunc = @(net, X, targetData, filterFunc)NNFilterLoss(net, X, targetData, filterFunc, false);
gradFunc = @(net, X, targetData, filterFunc)CalculateSalientMapping(net, X, targetData, filterFunc);

inputNames = {'wI', 'bA', 'mA', 'wL', 'rS', 'rR', 'zS', 'zR'};
for i = 1:numNetworks
    if matches(networks.filterType{i}, 'IIR')
        filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
        numOutputs = 2 * numFilters + 1;
    else
        filterFunc = @(output, target) BiquadLoss(output, target, numFilters, nfft, fs, fidx);
        numOutputs = 4 * numFilters + 1;
    end
    load([loadDir, filesep, listing{i}], 'loss', 'net', 'epochLosses', 'losses');
    figTitle = [networks.filterType{i}, ' Size: ', num2str(networks.networkSize{i}), ' Layers: ', num2str(networks.numLayers{i})];

    % PlotNNTrainingLossess(losses, epochLosses, figTitle)

    X = dlarray(single(inputData), "CB");
    gradients = dlfeval(gradFunc, net, X, targetData, filterFunc);
    [loss, ~, ~, prediction] = dlfeval(lossFunc, net, X, targetData, filterFunc);
    iLosses = sum(abs(prediction - targetData), 1) / numFreq;
    iSqLosses = sum((prediction - targetData) .^ 2, 1) / numFreq;

    mid = median(iLosses);
    av = mean(iLosses);
    percentileBounds = 5:5:95;
    percentiles = prctile(iLosses, percentileBounds);
    
    examplesIdxLo = find(iLosses < percentiles(1));
    examplesIdxHi = find(iLosses > percentiles(end));
    [~, examplesIdxMax] = maxk(iLosses, numSelect);
    
    idx = randi(length(examplesIdxLo), [1, numSelect]);
    examplesIdxLo = examplesIdxLo(idx);

    idx = randi(length(examplesIdxHi), [1, numSelect]);
    examplesIdxHi = examplesIdxHi(idx);

%     figure
%     colororder(colours)
%     semilogx(fc, targetData(:, examplesIdxLo))
%     hold on
%     semilogx(fc, prediction(:, examplesIdxLo), '--')
%     title(['below ', num2str(percentiles(1)), ' error: ', figTitle])
%     legend(num2str(inputData(:, examplesIdxLo)'), 'Location', 'southoutside')
%     xlim([20 20e3])
%     ylim([-60 20])

%     figure
%     colororder(colours)
%     semilogx(fc, targetData(:, examplesIdxHi))
%     hold on
%     semilogx(fc, prediction(:, examplesIdxHi), '--')
%     title(['above ', num2str(percentiles(end)), ' error: ', figTitle])
%     legend(num2str(inputData(:, examplesIdxHi)'), 'Location', 'southoutside')
%     xlim([20 20e3])
%     ylim([-60 20])

%     figure
%     colororder(colours)
%     semilogx(fc, targetData(:, examplesIdxMax))
%     hold on
%     semilogx(fc, prediction(:, examplesIdxMax), '--')
%     title(['Worst cases: ', figTitle])
%     legend(num2str(inputData(:, examplesIdxMax)'), 'Location', 'southoutside')
%     xlim([20 20e3])
%     ylim([-60 20])

%     figure
%     histogram(iLosses)
%     title(['Losses: ', figTitle])

    n = 3;
    m = 4;
    x = [];
    for k = 1:n * m
        for j = 1:length(gradients)
            x{k}(:,j) = gradients{j}(:,k);
        end
    end

    figure
    ax = gca;
    t = tiledlayout(n,m,'TileSpacing','Compact');
    for k = 1:n * m
        nexttile
        image(extractdata(x{k}), 'CDataMapping','scaled')
        clim([-1 1])
        colorbar
        for j = 1:8
            labels{j} = [inputNames{j}, ': ', num2str(inputData(j, k))];
        end
        yticks(1:8)
        yticklabels(labels)
        title('Loss: ', num2str(iLosses(k)))
    end
    title(t,[figTitle, ' Loss: ', num2str(percentiles(end))])
    xlabel(t,'Outputs')
    ylabel(t,'Inputs')
    
    allPercentiles(i,:) = percentiles;
%     figure
%     histogram(iSqLosses)
%     title(['Square losses: ', figTitle])
end



%% End

figure
plot(percentileBounds, allPercentiles)
ylim([0 6])

% cost = [];
% runningMean = [];
% for i = 1:numEpochs
%     idx = max(i - 5, 1):i - 1;
%     runningMean(i) = mean(epochLosses(idx));
%     cost(i) = runningMean(i) - epochLosses(i);
% end
% 
% xEpoch = 1:numEpochs;
% figure
% plot(xEpoch, runningMean)
% 
% figure
% plot(xEpoch, cost)

%% Best

close all

[~, i] = min(allPercentiles(:, end));

if matches(networks.filterType{i}, 'IIR')
    filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
    numOutputs = 2 * numFilters + 1;
else
    filterFunc = @(output, target) BiquadLoss(output, target, numFilters, nfft, fs, fidx);
    numOutputs = 4 * numFilters + 1;
end
load([loadDir, filesep, listing{i}], 'loss', 'net', 'epochLosses', 'losses');
figTitle = [networks.filterType{i}, ' Size: ', num2str(networks.networkSize{i}), ' Layers: ', num2str(networks.numLayers{i})];

% PlotNNTrainingLossess(losses, epochLosses, figTitle)

isWeight = matches(net.Learnables.Parameter, "Weights");

weights = net.Learnables.Value(isWeight);

numLayers = length(weights);

allWeights = [];
for i = 1:numLayers
    allWeights = [allWeights; reshape(extractdata(weights{i}), [], 1)];
end
threshold = 0;
weightThreshold = prctile(abs(allWeights), threshold);

for i = 1:numLayers
    weights{i}(abs(weights{i}) < weightThreshold) = 0;
end

net.Learnables.Value(isWeight) = weights;

X = dlarray(single(inputData), "CB");
gradients = dlfeval(gradFunc, net, X, targetData, filterFunc);
[loss, ~, ~, prediction] = dlfeval(lossFunc, net, X, targetData, filterFunc);
iLosses = sum(abs(prediction - targetData), 1) / numFreq;
iSqLosses = sum((prediction - targetData) .^ 2, 1) / numFreq;

% analyzeNetwork(net,X)
layers = {'InputLayer', 'ActivationLayer_1', 'ActivationLayer_2', 'ActivationLayer_3', 'ActivationLayer_4', 'ActivationLayer_5', 'OutputLayer'};
[inOutput, output1, output2, output3, output4, output5, outOutput] = predict(net, X, 'Outputs', layers);

figure
ax = gca;
t = tiledlayout(n,m,'TileSpacing','Compact');
image(extractdata(inOutput), 'CDataMapping','scaled')
%clim([-1 1])
colorbar
title(t,'Input Layer')
xlabel(t,'Outputs')
ylabel(t,'Inputs')

figure
ax = gca;
t = tiledlayout(n,m,'TileSpacing','Compact');
image(extractdata(output1), 'CDataMapping','scaled')
%clim([-1 1])
colorbar
title(t,'Layer 1')
xlabel(t,'Outputs')
ylabel(t,'Inputs')

figure
ax = gca;
t = tiledlayout(n,m,'TileSpacing','Compact');
image(extractdata(output2), 'CDataMapping','scaled')
%clim([-1 1])
colorbar
title(t,'Layer 2')
xlabel(t,'Outputs')
ylabel(t,'Inputs')

figure
ax = gca;
t = tiledlayout(n,m,'TileSpacing','Compact');
image(extractdata(output3), 'CDataMapping','scaled')
%clim([-1 1])
colorbar
title(t,'Layer 3')
xlabel(t,'Outputs')
ylabel(t,'Inputs')

figure
ax = gca;
t = tiledlayout(n,m,'TileSpacing','Compact');
image(extractdata(output4), 'CDataMapping','scaled')
%clim([-1 1])
colorbar
title(t,'Layer 4')
xlabel(t,'Outputs')
ylabel(t,'Inputs')

figure
ax = gca;
t = tiledlayout(n,m,'TileSpacing','Compact');
image(extractdata(output5), 'CDataMapping','scaled')
%clim([-1 1])
colorbar
title(t,'Layer 5')
xlabel(t,'Outputs')
ylabel(t,'Inputs')

figure
ax = gca;
t = tiledlayout(n,m,'TileSpacing','Compact');
image(extractdata(outOutput), 'CDataMapping','scaled')
%clim([-1 1])
colorbar
title(t,'Output Layer')
xlabel(t,'Outputs')
ylabel(t,'Inputs')

%% Test

[z, p, k] = CreateIIRFromNNOutput(outOutput, 2);

[tfmag, ~] = CreateIIRFilter(z, p, k, nfft, fs);
    
tfmagNBand = CreateNBandMagnitude(tfmag, fidx);

targetTfmagNBand = max(-128, min(128, tfmagNBand));


beta0 = zeros(8, 17);
func = @(beta, x)test(beta, x);

options = statset('nlinfit');
options.MaxIter = 10000;
z = double(extractdata(z));
p = extractdata(p);
k = extractdata(k);
[betaz1,R1,J1,CovB1] = nlinfit(inputData,z(1,:),func, beta0, 'options', options);
[betaz2,R2,J2,CovB2] = nlinfit(inputData,z(2,:),func, beta0, 'options', options);
[betap1,R3,J3,CovB3] = nlinfit(inputData,p(1,:),func, beta0, 'options', options);
[betap2,R4,J4,CovB4] = nlinfit(inputData,p(2,:),func, beta0, 'options', options);
[betak,R5,J5,CovB5] = nlinfit(inputData,k,func, beta0, 'options', options);

zPred = [test(betaz1, inputData); test(betaz2, inputData)];
pPred = [test(betap1, inputData); test(betap2, inputData)];
kPred = test(betak, inputData);

[tfmag, ~] = CreateIIRFilter(zPred, pPred, kPred, nfft, fs);
    
tfmagNBand = CreateNBandMagnitude(tfmag, fidx);

tfmagNBand = max(-128, min(128, tfmagNBand));

lossAbs = sum((tfmagNBand - targetData).^2, 'all')  / numel(tfmagNBand);
lossNN = sum((tfmagNBand - targetTfmagNBand).^2, 'all')  / numel(tfmagNBand);

mid = median(iLosses);
av = mean(iLosses);
percentileBounds = 5:5:95;
percentiles = prctile(iLosses, percentileBounds);

examplesIdxLo = find(iLosses < percentiles(1));
examplesIdxHi = find(iLosses > percentiles(end));
[~, examplesIdxMax] = maxk(iLosses, numSelect);

idx = randi(length(examplesIdxLo), [1, numSelect]);
examplesIdxLo = examplesIdxLo(idx);

idx = randi(length(examplesIdxHi), [1, numSelect]);
examplesIdxHi = examplesIdxHi(idx);

figure
colororder(colours)
semilogx(fc, targetData(:, examplesIdxLo))
hold on
semilogx(fc, prediction(:, examplesIdxLo), '--')
title(['below ', num2str(percentiles(1)), ' error: ', figTitle])
legend(num2str(inputData(:, examplesIdxLo)'), 'Location', 'southoutside')
xlim([20 20e3])
ylim([-60 20])

figure
colororder(colours)
semilogx(fc, targetData(:, examplesIdxHi))
hold on
semilogx(fc, prediction(:, examplesIdxHi), '--')
title(['above ', num2str(percentiles(end)), ' error: ', figTitle])
legend(num2str(inputData(:, examplesIdxHi)'), 'Location', 'southoutside')
xlim([20 20e3])
ylim([-60 20])

figure
colororder(colours)
semilogx(fc, targetData(:, examplesIdxMax))
hold on
semilogx(fc, prediction(:, examplesIdxMax), '--')
title(['Worst cases: ', figTitle])
legend(num2str(inputData(:, examplesIdxMax)'), 'Location', 'southoutside')
xlim([20 20e3])
ylim([-60 20])

figure
histogram(iLosses)
title(['Losses: ', figTitle])

n = 3;
m = 4;
x = [];
for k = 1:n * m
    for j = 1:length(gradients)
        x{k}(:,j) = gradients{j}(:,k);
    end
end

figure
ax = gca;
t = tiledlayout(n,m,'TileSpacing','Compact');
for k = 1:n * m
    nexttile
    image(extractdata(x{k}), 'CDataMapping','scaled')
    clim([-1 1])
    colorbar
    for j = 1:8
        labels{j} = [inputNames{j}, ': ', num2str(inputData(j, k))];
    end
    yticks(1:8)
    yticklabels(labels)
    title('Loss: ', num2str(iLosses(k)))
end
title(t,[figTitle, ' Loss: ', num2str(percentiles(end))])
xlabel(t,'Outputs')
ylabel(t,'Inputs')
    
%     figure
%     histogram(iSqLosses)
%     title(['Square losses: ', figTitle])


%% Export neural network

CheckFileDir('NNonnxExport')
exportONNXNetwork(net, ['NNonnxExport', filesep, 'testExport.onnx'])

function Y = test(beta, x)

    [numInputs, numData] = size(x);
    squares = zeros(numInputs, numInputs, numData);
    for i = 1:numData
        squares(:, :, i) = x(:,i) * x(:,i)';
    end
    for i = 1:numData
        triples(:, :, i) = squares(:, :, i) .* x(:,i);
    end
    Y = sum(beta(:,1) .* x, 1) + squeeze(sum(sum(beta(:,2:numInputs + 1) .* squares, 1),2))' + squeeze(sum(sum(beta(:,numInputs + 2:end) .* triples, 1),2))';
end