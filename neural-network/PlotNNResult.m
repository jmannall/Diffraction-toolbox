
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
inputData = dlarray(single(inputData), "CB");
targetData = extractdata(targetData);

numFilters = 2;
nBands = 12;
[~, tfmag, ~, fvec, ~] = DefaultBTM(controlparameters);
[~, fc, fidx] = CreateFrequencyNBands(tfmag, fvec, nBands);
numFreq = fidx(end);

%% Plots

close all

n = 5;
colours = colororder;
colours = colours(1:n, :);

%for i = 1:numNetworks
%     if matches(networks.filterType, 'IIR')
%         filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
%         numOutputs = 2 * numFilters + 1;
%     else
%         filterFunc = @(output, target) BiquadLoss(output, target, numFilters, nfft, fs, fidx);
%         numOutputs = 4 * numFilters + 1;
%     end
%     load([loadDir, filesep, listing{i}], 'loss', 'net', 'epochLosses', 'losses');
%     figTitle = [networks.filterType{i}, ' Size: ', num2str(networks.networkSize{i}), ' Layers: ', num2str(networks.numLayers{i})];
    filterFunc = @(output, target) BiquadLoss(output, target, numFilters, nfft, fs, fidx);
    numOutputs = 4 * numFilters + 1;
%     filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
%     numOutputs = 2 * numFilters + 1;
    figTitle = 'Biquad_1';
    PlotNNTrainingLossess(losses, epochLosses, figTitle)
    [loss, ~, ~, prediction] = NNFilterLoss(net, inputData, targetData, filterFunc, false);
    iLosses = sum(abs(prediction - targetData), 1) / numFreq;
    iSqLosses = sum((prediction - targetData) .^ 2, 1) / numFreq;

    mid = median(iLosses);
    av = mean(iLosses);
    percentiles = prctile(iLosses, [5 10 25 75 90 95]);
    
    examplesIdxLo = find(iLosses < percentiles(1), n);
    examplesIdxHi = find(iLosses > percentiles(6), n);
    [~, examplesIdxMax] = maxk(iLosses, n);
    
    figure
    colororder(colours)
    semilogx(fc, targetData(:, examplesIdxLo))
    hold on
    semilogx(fc, prediction(:, examplesIdxLo), '--')
    title(['below ', num2str(percentiles(1)), ' error: ', figTitle])
    legend(num2str(extractdata(inputData(:, examplesIdxLo))'), 'Location', 'southoutside')
    xlim([20 20e3])
    ylim([-60 20])

    figure
    colororder(colours)
    semilogx(fc, targetData(:, examplesIdxHi))
    hold on
    semilogx(fc, prediction(:, examplesIdxHi), '--')
    title(['above ', num2str(percentiles(6)), ' error: ', figTitle])
    legend(num2str(extractdata(inputData(:, examplesIdxHi))'), 'Location', 'southoutside')
    xlim([20 20e3])
    ylim([-60 20])

    figure
    colororder(colours)
    semilogx(fc, targetData(:, examplesIdxMax))
    hold on
    semilogx(fc, prediction(:, examplesIdxMax), '--')
    title(['Worst cases: ', figTitle])
    legend(num2str(extractdata(inputData(:, examplesIdxMax))'), 'Location', 'southoutside')
    xlim([20 20e3])
    ylim([-60 20])

    figure
    histogram(iLosses)
    title(['Losses: ', figTitle])

%     figure
%     histogram(iSqLosses)
%     title(['Square losses: ', figTitle])
%end


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
