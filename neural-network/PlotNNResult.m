
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

    figure
    histogram(iLosses)
    title(['Losses: ', figTitle])

    n = 3;
    m = 4;
    for k = 1:n * m
        for j = 1:length(gradients)
            x{k}(:,j) = gradients{j}(:,k);
        end
    end

    figure
    ax = gca;
    tiledlayout(n, m)
    for k = 1:n * m
        nexttile
        image(extractdata(x{k}), 'CDataMapping','scaled')
        clim([-1 1])
        colorbar
        xlabel('Outputs')
        ylabel('Inputs')
    end
    
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

%% Export neural network

CheckFileDir('NNonnxExport')
exportONNXNetwork(net, ['NNonnxExport', filesep, 'testExport.onnx'])
