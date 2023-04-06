
fs = 48e3;
nfft = 8192;
c = 344;
controlparameters = struct('fs', 2 * fs, 'nfft', 2 * nfft, 'difforder', 1, 'c', c, 'saveFiles', 2, 'noDirect', true);

epochSize = 20e3;

numInvalid = 0;
numMissing = 0;
pattern = ' ';
missing = pattern;
invalid = pattern;
numSets = 500;
fileStem = ['results' filesep 'CreateBtmTrainingData_PhiiFlat' filesep 'CreateBtmTrainingData_'];
for i = 1:numSets
    fileExist = exist([fileStem, num2str(i), '-48_NN.mat'], 'file');
    if fileExist == 2
        load([fileStem, num2str(i), '-48_NN.mat'])
    else
        missing = append(missing, [num2str(i) ', ']);
        numMissing = numMissing + 1;
    end
    excess(i) = sum(targetData(1,:) > 10);
    if excess(i) > 0
        invalid = append(invalid, [num2str(i) ', ']);
        numInvalid = numInvalid + 1;
    end
end

idx = {'TestData', 'ValidationData'};
for i = 1:length(idx)
    fileExist = exist([fileStem, idx{i}, '-48_NN.mat'], 'file');
    
    if fileExist == 2
        load([fileStem, idx{i}, '-48_NN.mat'])
    else
        missing = append(missing, [idx{i} ', ']);
        numMissing = numMissing + 1;
    end
    excess(numSets + i) = sum(targetData(1,:) > 10);
    if excess(i) > 0
        invalid = append(invalid, [idx{i} ', ']);
        numInvalid = numInvalid + 1;
    end
end
totalInvalid = sum(excess);

%% Result

if missing == pattern
    if invalid == pattern
        disp('Training data set is complete and valid')
    else
        disp(['Invalid:', invalid])
    end
else
    disp(['Missing:', missing])
    if invalid == pattern
        disp('No invalid data')
    else
        disp('No missing data')
    end
end
