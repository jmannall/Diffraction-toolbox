close all
clear all

directory = 'tempNN';
saveDir = 'tempNNNew';
CheckFileDir(saveDir)
dirInfo = dir(directory);
dirInfo = dirInfo(3:end);
disp('Start')

numFiles = length(dirInfo);
for j = 1:numFiles
    loadPath = [directory filesep dirInfo(j).name];
    load(loadPath)
    disp(['Load path: ' loadPath])
    numInputs = nP.numInputs;
    numOutputs = nP.numOutputs;
    alpha = nP.alpha;
    savePath = nP.savePath;
    saveDir = cell2mat(extractBetween(savePath,'/','/'));

    seed = nP.seed;
    nP = struct('numInputs', numInputs, 'numOutputs', numOutputs, 'alpha', alpha, 'saveDir', saveDir, 'savePath', savePath, 'seed', seed);

    idx = CreateNNIdx(hP, tP, nP);
    savePath = [saveDir filesep idx '.mat'];
    disp(['Save path: ' savePath])
    save(savePath, "net", "losses", "hP", "tP", "nP", "iteration", "i", '-v7.3')
    delete(loadPath)
end
disp('Complete')