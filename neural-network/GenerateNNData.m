function GenerateNNData(idx, fs)

    disp(idx)
    disp('Hello World')

    mkdir('results')
    mkdir(['results', filesep, 'CreateBtmTrainingData'])
    mkdir(['results', filesep, 'SingleWedge'])
    mkdir(['results', filesep, 'DefaultBTM'])

    %poolobj = gcp
    %if isempty(poolobj)
    %    parpool([1 48]);
    %end
    parpool([1 48]);

    epochSize = 20e3;
    
    nfft = 16384;
    c = 344;
    controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2, 'noDirect', true);
    
    % gpuDevice(1)
    seed = idx;
    rng(seed)
    
    seed = round(2 ^ 32 * rand(1)) - 1
    rng(seed)
    
    numDataSets = 10;

    for i = 1:numDataSets
        saveIdx = idx * numDataSets + i;
        [~, ~] = CreateBtmTrainingData(epochSize, controlparameters, saveIdx);
    end
end