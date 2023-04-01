
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
for i = 1:numNets
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

[lossA, idxA] = min(loss);
[bestLoss, idxB] = min(lossA);

disp(['The best network is ', char(networkNames(idxA(idxB))), ', Run ', num2str(idxB)])
