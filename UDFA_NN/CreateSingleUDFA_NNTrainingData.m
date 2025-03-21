function [trainingData, targetData, validationData] = CreateSingleUDFA_NNTrainingData(numInputs, controlparameters, saveValidation, index)

    fc = RandomLoguniformDistribution([1 10^7], numInputs);

    fs = controlparameters.fs;
    nfft = controlparameters.nfft;
    nBands = controlparameters.nBands;
    fvec = controlparameters.fvec;

    config = DataHash({numInputs, fs, nfft, nBands});

    % Save paths
    rootDir = 'NNData';
    saveDir = 'SingleUDFA_NN';
    CheckFileDir(rootDir)
    CheckFileDir([rootDir filesep saveDir])
    % Create file info
    if nargin < 3
        saveFile = ['Data-', num2str(config)];
    else
        saveFile = ['Data-', num2str(index), '-', num2str(config)];
    end
    savePath = [rootDir filesep saveDir filesep saveFile];

    complete = exist([cd filesep savePath '.mat'], "file");

    if complete == 2
        if (saveValidation)
            load([cd filesep savePath '.mat'], 'trainingData', 'targetData', 'validationData');
        else
            load([cd filesep savePath '.mat'], 'trainingData', 'targetData');
        end
    else
        disp('Generate data')
        % Generate data
        validationData = zeros(nfft / 2, numInputs);
        for i = 1:numInputs
            response = filterApprox(fvec, fc(i), 1, 0, 1.44, 0.2, 0.5);
            validationData(:,i) = mag2db(abs(response));
        end
        targetData = CreateFrequencyNBands(validationData, fvec, nBands);
        targetData = max(-128, min(128, targetData));
        trainingData = log10(fc)';
        if (saveValidation)
            save([cd filesep savePath '.mat'], 'trainingData', 'targetData', 'validationData');
        else
            save([cd filesep savePath '.mat'], 'trainingData', 'targetData');
        end
    end
end