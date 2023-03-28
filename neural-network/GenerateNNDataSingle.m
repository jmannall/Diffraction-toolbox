function GenerateNNDataSingle(size, fs)

    rng shuffle
    rng

    disp(size)
    disp('Hello World')

    mkdir('results')
    mkdir(['results', filesep, 'CreateBtmTrainingData'])
    mkdir(['results', filesep, 'SingleWedge'])
    mkdir(['results', filesep, 'DefaultBTM'])
    
    nfft = 16384;
    c = 344;
    controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);
    
    epochSize = 20e3;
    [~, ~] = CreateBtmTrainingData(epochSize, controlparameters);
end