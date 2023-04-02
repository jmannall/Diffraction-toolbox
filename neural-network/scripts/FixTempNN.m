directory = 'tempNN';
dirInfo = dir(directory);
dirInfo = dirInfo(3:end);
disp('Start')

numFiles = length(dirInfo);
for j = 1:numFiles
    loadPath = [directory filesep dirInfo(j).name];
    load(loadPath)

    idx = CreateNNIdx(hP, tP, nP);
    savePath = [directory filesep idx '.mat'];
    save(savePath, "net", "losses", "hP", "tP", "nP", "iteration", "i", '-v7.3')
    delete(loadPath)
end
disp('Complete')