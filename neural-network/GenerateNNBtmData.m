
close all

numEpochs = 500;
epochSize = 20e3;

fs = 96e3;
nfft = 16384;
c = 344;
controlparameters = struct('fs', fs, 'nfft', nfft, 'difforder', 1, 'c', c, 'saveFiles', 2);

gpuDevice(1)

for i = 1:numEpochs
    [trainingData, targetData] = CreateBtmTrainingData(epochSize, controlparameters, i);
end