
close all
clear all

colorStore = colororder;

%% Load networks

loadDir = 'NNSaves/SingleUDFA_NN';
dirInfo = dir(loadDir);
dirInfo = dirInfo(3:end);
idx = [dirInfo.isdir];
folders = dirInfo([dirInfo.isdir]);

[~, runIdx] = sort(str2double(extractAfter({folders.name}, 'Run')));
folders = folders(runIdx);
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

plotGraphs = false;
numNets = length(nets);
gx = 5;

color = colorStore(1:2,:);
for i = 1:numNets
    close all
    numRuns = length(nets{i});
    for j = 1:numRuns
        if (isempty(nets{i}(j).losses))
            continue;
        end
        numLayers = nets{i}(j).hP.numLayers;
        hiddenLayerSize = nets{i}(j).hP.hiddenLayerSize;
        learnRate = nets{i}(j).hP.learnRate;
        size = CalculateNNIIRCost(numLayers, hiddenLayerSize, nets{i}(j).nP.numInputs, nets{i}(j).nP.numOutputs, gx);
        if plotGraphs
            titleText = ['Net: ', num2str(numLayers), '-', num2str(hiddenLayerSize), ', Learn rate: ', num2str(learnRate), ', Size: ', num2str(size)];
            figure
            hold on
            colororder(color)
        end
        losses = nets{i}(j).losses;
        numEpochs = nets{i}(j).tP.numEpochs;
        lvec = 1:numEpochs;
        if plotGraphs
            plot(lvec, losses.test)
            plot(lvec, losses.epoch, '--')
        end
        numIterations = length(losses.iteration);
        %lvec = (1:numIterations) / numIterations * numEpochs;
        %plot(lvec, losses.iteration)
        sizes(i,j) = size;
        loss(i,j) = losses.test(end);
        % networkNames{i,j} = ['Net: ', num2str(numLayers), '-', num2str(hiddenLayerSize), ', Learn rate: ', num2str(learnRate), ', Size: ', num2str(size), ', Seed: ', num2str(nets{i}(j).nP.seed), ', Run: ', num2str(j), ', Loss: ', num2str(losses.test(end))];
        networkNames{i,j} = ['Net: ', num2str(numLayers), '-', num2str(hiddenLayerSize), ', Size: ', num2str(size), ', Seed: ', num2str(nets{i}(j).nP.seed), ', Run: ', num2str(j), ', Loss: ', num2str(losses.test(end))];
    end
    if plotGraphs
        grid on
        ylim([0 100])      
        title(titleText)
        ylabel('Epoch')
        xlabel('Mean absolute error (dB)')
    end
end

%%

[lossSort, b] = sort(reshape(loss, [], 1));
networkNamesReshape = reshape(networkNames, [], 1);
sizesReshape = reshape(sizes, [], 1);

%b = mod(b - 1, i) + 1;
networkNamesSort = networkNamesReshape(b);
sizesSort = sizesReshape(b);
% loss(loss == 0) = 1000;
% [lossA, idxA] = min(loss);
% [bestLoss, idxB] = min(lossA);

disp(['The best network is ', char(networkNamesSort(1))])

%%
close all

f = figure;
scatter(lossSort, sizesSort)
grid on
xlabel('Accuracy')
ylabel('Network Size')
title('Model size vs Accuracy')

dcm = datacursormode(f);
dcm.UpdateFcn = {@customDataTip, networkNamesSort};
function txt = customDataTip(~, eventData, IDs)
    s = eventData.Target;
    pos = eventData.Position;
    ID = IDs{s.XData == pos(1) & s.YData == pos(2)};
    valueFormat = ' \color[rgb]{0 0.6 1}\bf';
    removeValueFormat = '\color[rgb]{.25 .25 .25}\rm';
    txt = {['X',[valueFormat num2str(pos(1)) removeValueFormat]],...
           ['Y',[valueFormat num2str(pos(2)) removeValueFormat]],...
           ['ID',[valueFormat num2str(ID) removeValueFormat]]};
end