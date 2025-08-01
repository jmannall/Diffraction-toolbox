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

    if (complete == 2)
        if (saveValidation)
            load([cd filesep savePath '.mat'], 'trainingData', 'targetData', 'validationData');
        else
            load([cd filesep savePath '.mat'], 'trainingData', 'targetData');
        end
        trainingData(1:numNNInputs / 2,:) = log10(trainingData(1:numNNInputs / 2,:));
    else
        disp('Generate data')
        % Generate data
        validationData = zeros(nfft / 2, numInputs);
        trainingData = zeros(numNNInputs, numInputs);
        %fc = zeros(5, numInputs);
        %gain = zeros(5, numInputs);
        %blendExpn = zeros(4, numInputs);
        %Q = zeros(4, numInputs);
        %udfaTerms = zeros(70, 5, numInputs);
        %fAxis = logspace(log10(20), log10(2e4), numNNInputs);
        for i = 1:length(geometry.wedgeIndex)
            thetaW = deg2rad(geometry.wedgeIndex(i));
            thetaS = deg2rad(geometry.thetaS(i));
            thetaR = deg2rad(geometry.thetaR(i));
            zW = geometry.wedgeLength(i);
            rS = geometry.rS(i);
            rR = geometry.rR(i);
            zS = geometry.zS(i);
            zR = geometry.zR(i);
            %offset = mag2db(sqrt((rR+rS)^2+(zR-zS)^2));
            % Test run using infinite edges (simpler) only needs 2 x (fc, gain)
            % Finite edges adds 4 x (fc, gain), variable blendExpn and Q
            % zA not on edge makes 2 x (fc, gain) inf and adds a lpf with (fc, gain, set blendExpn and Q)
            [pDiffr, ~, ~, ~, pDiffr_terms, fcOut, gainOut, blendOut] = getUDFA(fvec,thetaW,[rS zS],[rR zR],thetaS,thetaR,[0 zW],op);
            %geometry.wedgeLength(i) = inf;
            %[pDiffr, ~, ~, ~, pDiffr_terms, fcOut, gainOut, blendOut] = getUDFA(fvec,thetaW,[rS zS],[rR zR],thetaS,thetaR,[],op);
            % Remove distance factor
            distance = sqrt((rR+rS)^2+(zR-zS)^2);
            pDiffr = distance * pDiffr;
            validationData(:,i) = mag2db(abs(pDiffr));

            trainingData(:,i) = [fcOut(1:numNNInputs / 2), gainOut(1:numNNInputs / 2)];
            %fc(:,i) = fcOut';
            %gain(:,i) = gainOut';
            %blendExpn(:,i) = blendOut(1:4)';
            %Q(:,i) = blendOut(6:9)';

            % if (size(pDiffr_terms, 2) < 5)
            %     pDiffr_terms = [pDiffr_terms, zeros(nfft / 2, 1)];
            % end
            % udfaTerms(:,:,i) = CreateFrequencyNBands(pDiffr_terms, fvec, nBands);
    
            % OLD for if training NN using target freq response.
            % pDiffr = getUDFA(fAxis,thetaW,[rS zS],[rR zR],thetaS,thetaR,[0 zW],op);
            % trainingData(:,i) = mag2db(abs(pDiffr)) - offset;
        end
        targetData = CreateFrequencyNBands(validationData, fvec, nBands);
        targetData = max(-128, min(128, targetData));
        %trainingData = max(-128, min(128, trainingData)) / 128;
        if (saveValidation)
            save([cd filesep savePath '.mat'], 'trainingData', 'targetData', 'validationData');
        else
            save([cd filesep savePath '.mat'], 'trainingData', 'targetData');
        end
        %save([cd filesep savePath2 '.mat'], 'fc', 'gain', 'blendExpn', 'Q', 'udfaTerms');
    end
end