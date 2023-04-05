function GenerateNNDataSingle(size, fs)

    rng shuffle
    rng

    disp(size)
    disp('Hello World')

    CheckFileDir('results')
    CheckFileDir(['results', filesep, 'CreateBtmTrainingData'])
    CheckFileDir(['results', filesep, 'SingleWedge'])
    CheckFileDir(['results', filesep, 'DefaultBTM'])
    
    nfft = 16384;
    c = 344;
    controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2, 'noDirect', true);
    
    epochSize = 20e3;
    [~, ~] = CreateBtmTrainingData(epochSize, controlparameters);
end