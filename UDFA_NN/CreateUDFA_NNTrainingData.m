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
    saveDir2 = 'FREQ_NN';
    CheckFileDir(rootDir)
    CheckFileDir([rootDir filesep saveDir])
    CheckFileDir([rootDir filesep saveDir2])
    % Create file info
    if nargin < 3
        saveFile = ['Data-', num2str(config)];
    else
        saveFile = ['Data-', num2str(index), '-', num2str(config)];
    end
    savePath = [rootDir filesep saveDir filesep saveFile];
    savePath2 = [rootDir filesep saveDir2 filesep saveFile];

    complete = exist([cd filesep savePath '.mat'], "file");
    complete2 = exist([cd filesep savePath2 '.mat'], "file");

    if (complete == 2 && complete2 == 2)
        if (saveValidation)
            load([cd filesep savePath '.mat'], 'trainingData', 'targetData', 'validationData');
        else
            load([cd filesep savePath '.mat'], 'trainingData', 'targetData');
        end
    else
        disp('Generate data')
        % Generate data
        validationData = zeros(nfft / 2, numInputs);
        trainingData = zeros(numNNInputs, numInputs);
        fc = zeros(5, numInputs);
        gain = zeros(5, numInputs);
        blendExpn = zeros(4, numInputs);
        Q = zeros(4, numInputs);
        udfaTerms = zeros(70, 5, numInputs);
        fAxis = logspace(log10(20), log10(2e4), numNNInputs);
        for i = 1:length(geometry.wedgeIndex)
            thetaW = deg2rad(geometry.wedgeIndex(i));
            thetaS = deg2rad(geometry.thetaS(i));
            thetaR = deg2rad(geometry.thetaR(i) - 180);
            zW = geometry.wedgeLength(i);
            rS = geometry.rS(i);
            rR = geometry.rR(i);
            zS = geometry.zS(i);
            zR = geometry.zR(i);
            %offset = mag2db(sqrt((rR+rS)^2+(zR-zS)^2));
            [pDiffr, ~, ~, ~, pDiffr_terms, fcOut, gainOut, blendOut] = getUDFA(fvec,thetaW,[rS zS],[rR zR],thetaS,thetaR,[0 zW],op);
            offset = mag2db(abs(pDiffr(1)));
            validationData(:,i) = mag2db(abs(pDiffr)) - offset;

            fc(:,i) = fcOut';
            gain(:,i) = gainOut';
            blendExpn(:,i) = blendOut(1:4)';
            Q(:,i) = blendOut(6:9)';

            if (size(pDiffr_terms, 2) < 5)
                pDiffr_terms = [pDiffr_terms, zeros(nfft / 2, 1)];
            end
            udfaTerms(:,:,i) = CreateFrequencyNBands(pDiffr_terms, fvec, nBands);
    
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
        save([cd filesep savePath2 '.mat'], 'fc', 'gain', 'blendExpn', 'Q', 'udfaTerms');
    end
end