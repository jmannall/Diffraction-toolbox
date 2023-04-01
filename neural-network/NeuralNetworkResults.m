
close all
clear all

colorStore = colororder;

%% Load networks

loadDir = 'NNSaves';
dirInfo = dir(loadDir);
dirInfo = dirInfo(3:end);
idx = [dirInfo.isdir];
folders = dirInfo([dirInfo.isdir]);

numFolders = length(folders);

for i = 1:numFolders
    dirPath = [loadDir filesep folders(i).name];
    dirInfo = dir(dirPath);
    dirInfo = dirInfo(3:end);

    numNetworks = length(dirInfo);
    for j = 1:numNetworks
        filePath = [dirPath filesep dirInfo(j).name];
        load(filePath)
        nets{j}(i) = struct('net', net, 'losses', losses, 'hP', hP, 'tP', tP, 'nP', nP);
    end
end

%% Plot training loss

numNets = length(nets);
numRuns = numFolders;

color = colorStore(1:2,:);
for i = 11:numNets
    numLayers = nets{i}(1).hP.numLayers;
    hiddenLayerSize = nets{i}(1).hP.hiddenLayerSize;
    learnRate = nets{i}(1).hP.learnRate;
    titleText = ['Net: ', num2str(numLayers), '-', num2str(hiddenLayerSize), ', Learn rate: ', num2str(learnRate)];
    
    figure
    hold on
    colororder(color)
    for j = 1:numRuns
        losses = nets{i}(j).losses;
        numEpochs = nets{i}(j).tP.numEpochs;
        lvec = 1:numEpochs;
        plot(lvec, losses.test)
        plot(lvec, losses.epoch, '--')
        numIterations = length(losses.iteration);
        %lvec = (1:numIterations) / numIterations * numEpochs;
        %plot(lvec, losses.iteration)
        loss(i,j) = losses.test(end);
    end
    grid on
    ylim([0 200])      
    title(titleText)
    ylabel('Epoch')
    xlabel('Mean absolute error (dB)')

    networkNames(i,1) = string(titleText);
end
%%
nfft = 8192;
fs = 48e3;
c = 344;
controlparameters = struct('fs', 2 * fs, 'nfft', 2 * nfft, 'difforder', 1, 'c', c, 'saveFiles', 2, 'noDirect', true);
    
[input, target] = CreateBtmTrainingData(20e3, controlparameters, 'TestData');
input = dlarray(input, "CB");
% Filter parameters
numFilters = 2;
nBands = 8;
controlparameters.fs = fs;
controlparameters.nfft = nfft;
[~, tfmag, ~, fvec, ~] = DefaultBTM(controlparameters);
[~, ~, fidx] = CreateFrequencyNBands(tfmag, fvec, nBands);


filterFunc = @(output, target) IIRFilterLoss(output, target, numFilters, nfft, fs, fidx);
test = nets{1}.net;
[loss, state, gradients, predictionN, prediction] = NNFilterLoss(test, input, target, filterFunc, false);

%%

[lossA, idxA] = min(loss);
[bestLoss, idxB] = min(lossA);

disp(['The best network is ', char(networkNames(idxA(idxB))), ', Seed: ', num2str(nets{idxA(idxB)}(idxB).nP.seed), ', Run ', num2str(idxB)])