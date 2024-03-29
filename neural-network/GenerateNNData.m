function GenerateNNData(idx, fs)

    idxList = [1, 3, 35, 37, 39];
    idx = idxList(idx + 1);
    seed = idx;
    rng(seed)
    seed = round(2 ^ 32 * rand(1)) - 1
    rng(seed)

    disp(idx)
    disp('Hello World')

    CheckFileDir('results')
    CheckFileDir(['results', filesep, 'CreateBtmTrainingData'])
    CheckFileDir(['results', filesep, 'SingleWedge'])
    CheckFileDir(['results', filesep, 'DefaultBTM'])

    nfft = 16384;
    c = 344;
    controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2, 'noDirect', true);
    
    epochSize = 20e3;
    numDataSets = 10;
    for i = 1:numDataSets
        saveIdx = idx * numDataSets + i;
        [~, ~] = CreateBtmTrainingData(epochSize, controlparameters, saveIdx);
    end
end