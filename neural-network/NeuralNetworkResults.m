
%close all
%clear all
      
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
gx = 5;

color = colorStore(1:2,:);
for i = 1:numNets
    numLayers = nets{i}(1).hP.numLayers;
    hiddenLayerSize = nets{i}(1).hP.hiddenLayerSize;
    learnRate = nets{i}(1).hP.learnRate;
    size = CalculateNNIIRCost(numLayers, hiddenLayerSize, nets{i}(1).nP.numInputs, nets{i}(1).nP.numOutputs, gx);
    titleText = ['Net: ', num2str(numLayers), '-', num2str(hiddenLayerSize), ', Learn rate: ', num2str(learnRate), ', Size: ', num2str(size)];
    
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
    % Add network file path
end

[lossSort, b] = sort(reshape(loss, [], 1));

b = mod(b - 1, i) + 1;
networkNamesSort = networkNames(b);

loss(loss == 0) = 1000;
[lossA, idxA] = min(loss);
[bestLoss, idxB] = min(lossA);

disp(['The best network is ', char(networkNames(idxA(idxB))), ', Seed: ', num2str(nets{idxA(idxB)}(idxB).nP.seed), ', Run ', num2str(idxB), ', Loss: ', num2str(bestLoss)])

%%

cost = CalculateNNIIRCost(6, 24, 8, 5, 5);