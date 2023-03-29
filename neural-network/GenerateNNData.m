function GenerateNNData(idx, fs)

    seed = idx;
    rng(seed)
    seed = round(2 ^ 32 * rand(1)) - 1
    rng(seed)

    disp(idx)
    disp('Hello World')

    mkdir('results')
    mkdir(['results', filesep, 'CreateBtmTrainingData'])
    mkdir(['results', filesep, 'SingleWedge'])
    mkdir(['results', filesep, 'DefaultBTM'])

    nfft = 16384;
    c = 344;
    controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2, 'noDirect', true);
    
    epochSize = 20e3;
    numDataSets = 5;
    for i = 1:numDataSets
        saveIdx = idx * numDataSets + i;
        [~, ~] = CreateBtmTrainingData(epochSize, controlparameters, saveIdx);
    end
end