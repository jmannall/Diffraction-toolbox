
close all
clear all

%% Load networks

loadDir = 'NNSaves';
listingStore = dir(loadDir);
listingStore = listingStore(3:end);
store = {listingStore.name};
filterType = extractBefore(store, '-');
idx = [];
for i = 1:length(listingStore)
    idx(i) = matches(filterType{i}, "5_iir") || matches(filterType{i}, "5_iirW");
end
listingStore = listingStore(idx == 1);
listing = {listingStore.name};
networks.filterType = extractBefore(listing, '-');
networks.networkSize = extractBetween(listing, '-', '_');
networkStructure = extractBetween(listing, '099-', '.');
networks.numLayers = extractBefore(networkStructure, '-');
networks.numNodes = extractAfter(networkStructure, '-');

numNetworks = length(listing);
[numEpochs, allLoss] = deal(zeros(1, numNetworks));
for i = 1:numNetworks
    load([loadDir, filesep, listing{i}], 'loss', 'net', 'epochLosses', 'losses');
    numEpochs(i) = length(epochLosses);
    allLoss(i) = loss;
    firstLoss(i) = epochLosses(1);
    % PlotNNTrainingLossess(losses, epochLosses, listing{i})
end

maxFirstLoss = max(firstLoss);
 
%% Go
finishedTraining = numEpochs == 500;
wellTrained = allLoss < 20;
numTrained = sum(finishedTraining);
numWellTrained = sum(wellTrained);
idx = finishedTraining & wellTrained;
num = sum(idx);
listing = {listingStore(idx).name};

idx = finishedTraining;
listing2 = {listingStore(idx).name};

networkStructure = extractBetween(listing, '099-', '.');
layers1 = cell2mat(extractBefore(networkStructure, '-'));
networkStructure = extractBetween(listing2, '099-', '.');
layers2 = cell2mat(extractBefore(networkStructure, '-'));

size2 = extractBetween(listing2, '-', '_');
x = 1;
%% 

% close all
% 
% layers = 2:7;
% sizes = 2000:2000:10000;
% 
% loadDir = 'NNSaves';
% listingStore = dir(loadDir);
% listingStore = listingStore(3:end);
% store = {listingStore.name};
% filterType = extractBefore(store, '-');
% idx = [];
% for i = 1:length(listingStore)
%     idx(i) = matches(filterType{i}, "iir");
% end
% listingStore = listingStore(idx == 1);
% 
% [lossGrid, trainedGrid] = deal(zeros(length(layers), length(sizes)));
% for i = 1:numNetworks / 4
%     load([loadDir, filesep, listingStore(i).name], 'loss', 'net', 'epochLosses', 'losses');
%     layer = str2double(networks.numLayers{i});
%     size = round(str2double(networks.networkSize{i}) / 1e3) * 1e3;
%     lossGrid(layer - 1, size / 2e3) = loss;
%     numEpochs = length(epochLosses);
%     trainedGrid(layer - 1, size / 2e3) = numEpochs == 500;
% end
% 
% figure
% surf(trainedGrid, 'FaceColor','none','EdgeColor','interp')
% title('Finished training iir')
% xlabel('network size')
% ylabel('number of layers')
% xticks(1:5)
% yticks(1:6)
% xlim([1 5])
% ylim([1 6])
% xticklabels(sizes)
% yticklabels(layers)
% 
% figure
% surf(lossGrid, 'FaceColor','none','EdgeColor','interp')
% zlim([0 10])
% clim([0 10])
% title('End loss iir')
% xlabel('network size')
% ylabel('number of layers')
% xticks(1:5)
% yticks(1:6)
% xlim([1 5])
% ylim([1 6])
% xticklabels(sizes)
% yticklabels(layers)
% 
% loadDir = 'NNSaves';
% listingStore = dir(loadDir);
% listingStore = listingStore(3:end);
% store = {listingStore.name};
% filterType = extractBefore(store, '-');
% idx = [];
% for i = 1:length(listingStore)
%     idx(i) = matches(filterType{i}, "iirW");
% end
% listingStore = listingStore(idx == 1);
% 
% [lossGrid, trainedGrid] = deal(zeros(length(layers), length(sizes)));
% for i = 1:numNetworks / 4
%     load([loadDir, filesep, listingStore(i).name], 'loss', 'net', 'epochLosses', 'losses');
%     layer = str2double(networks.numLayers{i});
%     size = round(str2double(networks.networkSize{i}) / 1e3) * 1e3;
%     lossGrid(layer - 1, size / 2e3) = loss;
%     numEpochs = length(epochLosses);
%     trainedGrid(layer - 1, size / 2e3) = numEpochs == 500;
% end
% 
% figure
% surf(trainedGrid, 'FaceColor','none','EdgeColor','interp')
% title('Finished training iirW')
% xlabel('network size')
% ylabel('number of layers')
% xticks(1:5)
% yticks(1:6)
% xlim([1 5])
% ylim([1 6])
% xticklabels(sizes)
% yticklabels(layers)
% 
% figure
% surf(lossGrid, 'FaceColor','none','EdgeColor','interp')
% zlim([0 10])
% clim([0 10])
% title('End loss iirW')
% xlabel('network size')
% ylabel('number of layers')
% xticks(1:5)
% yticks(1:6)
% xlim([1 5])
% ylim([1 6])
% xticklabels(sizes)
% yticklabels(layers)
% 
% loadDir = 'NNSaves';
% listingStore = dir(loadDir);
% listingStore = listingStore(3:end);
% store = {listingStore.name};
% filterType = extractBefore(store, '-');
% idx = [];
% for i = 1:length(listingStore)
%     idx(i) = matches(filterType{i}, "2_iir");
% end
% listingStore = listingStore(idx == 1);
% 
% [lossGrid, trainedGrid] = deal(zeros(length(layers), length(sizes)));
% for i = 1:numNetworks / 4
%     load([loadDir, filesep, listingStore(i).name], 'loss', 'net', 'epochLosses', 'losses');
%     layer = str2double(networks.numLayers{i});
%     size = round(str2double(networks.networkSize{i}) / 1e3) * 1e3;
%     lossGrid(layer - 1, size / 2e3) = loss;
%     numEpochs = length(epochLosses);
%     trainedGrid(layer - 1, size / 2e3) = numEpochs == 500;
% end
% 
% figure
% surf(trainedGrid, 'FaceColor','none','EdgeColor','interp')
% title('Finished training 2 iir')
% xlabel('network size')
% ylabel('number of layers')
% xticks(1:5)
% yticks(1:6)
% xlim([1 5])
% ylim([1 6])
% xticklabels(sizes)
% yticklabels(layers)
% 
% figure
% surf(lossGrid, 'FaceColor','none','EdgeColor','interp')
% zlim([0 10])
% clim([0 10])
% title('End loss 2 iir')
% xlabel('network size')
% ylabel('number of layers')
% xticks(1:5)
% yticks(1:6)
% xlim([1 5])
% ylim([1 6])
% xticklabels(sizes)
% yticklabels(layers)
% 
% loadDir = 'NNSaves';
% listingStore = dir(loadDir);
% listingStore = listingStore(3:end);
% store = {listingStore.name};
% filterType = extractBefore(store, '-');
% idx = [];
% for i = 1:length(listingStore)
%     idx(i) = matches(filterType{i}, "2_iirW");
% end
% listingStore = listingStore(idx == 1);
% 
% [lossGrid, trainedGrid] = deal(zeros(length(layers), length(sizes)));
% for i = 1:numNetworks / 4
%     load([loadDir, filesep, listingStore(i).name], 'loss', 'net', 'epochLosses', 'losses');
%     layer = str2double(networks.numLayers{i});
%     size = round(str2double(networks.networkSize{i}) / 1e3) * 1e3;
%     lossGrid(layer - 1, size / 2e3) = loss;
%     numEpochs = length(epochLosses);
%     trainedGrid(layer - 1, size / 2e3) = numEpochs == 500;
% end
% 
% figure
% surf(trainedGrid, 'FaceColor','none','EdgeColor','interp')
% title('Finished training 2 iirW')
% xlabel('network size')
% ylabel('number of layers')
% xticks(1:5)
% yticks(1:6)
% xlim([1 5])
% ylim([1 6])
% xticklabels(sizes)
% yticklabels(layers)
% 
% figure
% surf(lossGrid, 'FaceColor','none','EdgeColor','interp')
% zlim([0 10])
% clim([0 10])
% title('End loss 2 iirW')
% xlabel('network size')
% ylabel('number of layers')
% xticks(1:5)
% yticks(1:6)
% xlim([1 5])
% ylim([1 6])
% xticklabels(sizes)
% yticklabels(layers)

%% Create Data

fs = 48e3;
nfft = 8192;
c = 344;
testSize = 20e3;

controlparameters = struct('fs', fs, 'nfft',nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);
numFilters = 2;
nBands = 8;
[~, tfmag, ~, fvec, ~] = DefaultBTM(controlparameters);
[~, fc, fidx] = CreateFrequencyNBands(tfmag, fvec, nBands);
numFreq = fidx(end);

controlparameters.fs = 2 * fs;
controlparameters.nfft = 2 * nfft;
[inputData, targetData] = CreateBtmTrainingData(testSize, controlparameters, testSize);

%% Plots

close all

numSelect = 5;
colours = colororder;
colours = colours(1:numSelect, :);
lossFunc = @(net, X, targetData, filterFunc)NNFilterLoss(net, X, targetData, filterFunc, false);
gradFunc = @(net, X, targetData, filterFunc)CalculateSalientMapping(net, X, targetData, filterFunc);

inputNames = {'wI', 'bA', 'mA', 'wL', 'rS', 'rR', 'zS', 'zR'};
for i = 1:num
    filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
    numOutputs = 2 * numFilters + 1;
%     if matches(networks.filterType{i}, '2_iir')
%         filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
%         numOutputs = 2 * numFilters + 1;
%     else
%         filterFunc = @(output, target) BiquadLoss(output, target, numFilters, nfft, fs, fidx);
%         numOutputs = 4 * numFilters + 1;
%     end
    load([loadDir, filesep, listing{i}], 'loss', 'net', 'epochLosses', 'losses');
    figTitle = [networks.filterType{i}, ' Size: ', num2str(networks.networkSize{i}), ' Layers: ', num2str(networks.numLayers{i})];

    PlotNNTrainingLossess(losses, epochLosses, figTitle)

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

%     n = 3;
%     m = 4;
%     x = [];
%     for k = 1:n * m
%         for j = 1:length(gradients)
%             x{k}(:,j) = gradients{j}(:,k);
%         end
%     end

%     figure
%     ax = gca;
%     t = tiledlayout(n,m,'TileSpacing','Compact');
%     for k = 1:n * m
%         nexttile
%         image(extractdata(x{k}), 'CDataMapping','scaled')
%         clim([-1 1])
%         colorbar
%         for j = 1:8
%             labels{j} = [inputNames{j}, ': ', num2str(inputData(j, k))];
%         end
%         yticks(1:8)
%         yticklabels(labels)
%         title('Loss: ', num2str(iLosses(k)))
%     end
%     title(t,[figTitle, ' Loss: ', num2str(percentiles(end))])
%     xlabel(t,'Outputs')
%     ylabel(t,'Inputs')
    
    allPercentiles(i,:) = percentiles;
%     figure
%     histogram(iSqLosses)
%     title(['Square losses: ', figTitle])
    disp(i)
end

[allPercentiles, order] = sortrows(allPercentiles, 19);
idx = allPercentiles(:,end) < 20;

%% End

close all

figure
plot(percentileBounds, allPercentiles(1:7, :))
grid on
ylim([0 4])
ylabel('error (dB)')
xlabel('percentile')
legend(listing{order(idx)}, 'Location', 'eastoutside')
% cost = [];
% runningMean = [];
% for i = 1:numEpochs
%     idx = max(i - 5, 1):i - 1;
%     runningMean(i) = mean(epochLosses(idx));
%     cost(i) = runningMean(i) - epochLosses(i);
% end
% Create
% xEpoch = 1:numEpochs;
% figure
% plot(xEpoch, runningMean)
% 
% figure
% plot(xEpoch, cost)
edges = 0:0.1:4;

figure
t = tiledlayout(3, 3, 'TileSpacing', 'compact');
for i = 1:9
    nexttile
    histogram(allPercentiles(i,:), edges)
    title(listing{order(i)})
end

%% Best

[~, i] = min(allPercentiles(:, end));

if matches(networks.filterType{i}, 'iir')
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
curLimit = 0;
for k = 1:n * m
    for j = 1:length(gradients)
        x{k}(:,j) = gradients{j}(:,k);
        curLimit = max(curLimit, max(extractdata(abs(x{k})), [], "all"));
    end
end

limits = [-curLimit, curLimit];
figure
ax = gca;
t = tiledlayout(n,m,'TileSpacing','Compact');
for k = 1:n * m
    nexttile
    image(extractdata(x{k}), 'CDataMapping','scaled')
    clim(limits)
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
