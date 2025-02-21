function [trainingData, targetData, validationData] = CreateUDFA_NNTrainingData(numInputs, controlparameters, saveValidation, index)

    geometry = RandomGeometryUDFA(numInputs);

    fs = controlparameters.fs;
    nfft = controlparameters.nfft;
    nBands = controlparameters.nBands;
    numNNInputs = controlparameters.numNNInputs;
    fvec = controlparameters.fvec;
    
    op.udfa_underlying = 'pierce';
    op.UDFAdoAsymBlendF = 0;

    config = DataHash({numInputs, fs, nfft, nBands, numNNInputs});

    % Save paths
    rootDir = 'NNData';
    saveDir = 'UDFA_NN';
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
        disp('Load data')
        if (saveValidation)
            load([cd filesep savePath '.mat'], 'trainingData', 'targetData', 'validationData');
        else
            load([cd filesep savePath '.mat'], 'trainingData', 'targetData');
        end
    else
        % Generate data
        validationData = zeros(nfft / 2, numInputs);
        trainingData = zeros(numNNInputs, numInputs);
        fAxis = logspace(log10(20), log10(2e4), numNNInputs);
        for i = 1:length(geometry.wedgeIndex)
            thetaW = geometry.wedgeIndex(i);
            thetaS = geometry.thetaS(i);
            thetaR = geometry.thetaR(i);
            zW = geometry.wedgeLength(i);
            rS = geometry.rS(i);
            rR = geometry.rR(i);
            zS = geometry.zS(i);
            zR = geometry.zR(i);
            %offset = mag2db(sqrt((rR+rS)^2+(zR-zS)^2));
            pDiffr = getUDFA(fvec,thetaW,[rS zS],[rR zR],thetaS,thetaR,[0 zW],op);
            offset = mag2db(abs(pDiffr(1)));
            validationData(:,i) = mag2db(abs(pDiffr)) - offset;
    
            pDiffr = getUDFA(fAxis,thetaW,[rS zS],[rR zR],thetaS,thetaR,[0 zW],op);
            trainingData(:,i) = mag2db(abs(pDiffr)) - offset;
        end
        targetData = CreateFrequencyNBands(validationData, fvec, nBands);
        targetData = max(-128, min(128, targetData));
        trainingData = max(-128, min(128, trainingData)) / 128;
        if (saveValidation)
            save([cd filesep savePath '.mat'], 'trainingData', 'targetData', 'validationData');
        else
            save([cd filesep savePath '.mat'], 'trainingData', 'targetData');
        end
    end
end