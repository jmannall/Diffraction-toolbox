
close all

numEpochs = 350;
epochSize = 20e3;

fs = 96e3;
nfft = 16384;
c = 344;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);

% gpuDevice(1)
rng(6000)

for i = 289:numEpochs
    [trainingData, targetData] = CreateBtmTrainingData(epochSize, controlparameters, i);
end