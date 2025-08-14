function [trainingData, targetData, validationData, geometry] = CreateUDFA_NNTrainingData(numInputs, controlparameters, saveValidation, index)
    
    geometry = RandomGeometryUDFA(numInputs);

    % geometry.wedgeIndex = 270;
    % geometry.bendingAngle = 190;
    % geometry.minAngle = 5;
    % geometry.wedgeLength = 10;
    % geometry.rS = 2;
    % geometry.rR = 2;
    % geometry.zS = 5;
    % geometry.zR = 5;
    % geometry.zA = 5;
    % geometry.thetaS = geometry.minAngle;
    % geometry.thetaR = geometry.minAngle + geometry.bendingAngle;
    % numInputs = length(geometry.wedgeIndex);

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
    % saveDir2 = 'FREQ_NN';
    CheckFileDir(rootDir)
    CheckFileDir([rootDir filesep saveDir])
    % CheckFileDir([rootDir filesep saveDir2])
    % Create file info
    if nargin < 3
        saveFile = ['Data-', num2str(config)];
    else
        saveFile = ['Data-', num2str(index), '-', num2str(config)];
    end
    savePath = [rootDir filesep saveDir filesep saveFile];
    % savePath2 = [rootDir filesep saveDir2 filesep saveFile];

    complete = exist([cd filesep savePath '.mat'], "file");
    % complete2 = exist([cd filesep savePath2 '.mat'], "file");

    if (complete == 2)
        if (saveValidation)
            load([cd filesep savePath '.mat'], 'trainingData', 'targetData', 'validationData', 'geometry');
        else
            load([cd filesep savePath '.mat'], 'trainingData', 'targetData');
        end
    else
        tic
        disp('Generate data')
        % Generate data
        validationData = zeros(nfft / 2, numInputs);
        trainingData = zeros(numNNInputs, numInputs);
        % fc = zeros(5, numInputs);
        % gain = zeros(5, numInputs);
        % blendExpn = zeros(4, numInputs);
        % Q = zeros(4, numInputs);
        % udfaTerms = zeros(70, 5, numInputs);
        % fAxis = logspace(log10(20), log10(2e4), numNNInputs);

        if isempty(gcp('nocreate'))
            parpool('Threads', [1 128])
        end
        wedgeIndexAll = deg2rad(geometry.wedgeIndex);
        thetaSAll = deg2rad(geometry.thetaS);
        thetaRAll = deg2rad(geometry.thetaR);
        zWAll = geometry.wedgeLength;
        rSAll = geometry.rS;
        rRAll = geometry.rR;
        zSAll = geometry.zS;
        zRAll = geometry.zR;
        parfor i = 1:length(geometry.wedgeIndex)
            thetaW = wedgeIndexAll(i);
            thetaS = thetaSAll(i);
            thetaR = thetaRAll(i);
            zW = zWAll(i);
            rS = rSAll(i);
            rR = rRAll(i);
            zS = zSAll(i);
            zR = zRAll(i);
            %offset = mag2db(sqrt((rR+rS)^2+(zR-zS)^2));
            % Infinite edges 2 x (fc, gain)
            % Finite edges 4 x (fc, gain) + variable blendExpn, Q (fitted
            % and depend on fc and gain - maybe not needed for network)
            % zA not on edge makes 2 x (fc, gain) inf and adds a lpf with
            % (fc, gain, set blendExpn and Q)
            [~, ~, ~, ~, ~, fcOut, gainOut, ~] = getUDFA(fvec,thetaW,[rS zS],[rR zR],thetaS,thetaR,[0 zW],op);

            pDiffrBTMS = edBTMS(fvec,thetaW,[rS zS],[rR zR],thetaS,thetaR,[0 zW]);
            % geometry.wedgeLength(i) = inf;
            % [pDiffrInf, ~, ~, ~, pDiffr_termsInf, fcOutInf, gainOutInf, blendOutInf] = getUDFA(fvec,thetaW,[rS zS],[rR zR],thetaS,thetaR,[],op);

            % Remove distance factor
            distance = sqrt((rR+rS)^2+(zR-zS)^2);
            pDiffrBTMS = distance * pDiffrBTMS;
            validationData(:,i) = mag2db(abs(pDiffrBTMS));

            [fcOut, idx] = sort(fcOut(1:numNNInputs / 2));
            gainOut = gainOut(idx);

            trainingData(:,i) = [fcOut, gainOut];
            %fc(:,i) = fcOut';
            %gain(:,i) = gainOut';
            %blendExpn(:,i) = blendOut(1:4)';
            %Q(:,i) = blendOut(6:9)';

            % if (size(pDiffr_terms, 2) < 5)
            %     pDiffr_terms = [pDiffr_terms, zeros(nfft / 2, 1)];
            % end
            % udfaTerms(:,:,i) = CreateFrequencyNBands(pDiffr_terms, fvec, nBands);
    
            % OLD if training NN using target freq response
            % pDiffr = getUDFA(fAxis,thetaW,[rS zS],[rR zR],thetaS,thetaR,[0 zW],op);
            % offset = mag2db(abs(pDiffr(1)));
            % trainingData(:,i) = mag2db(abs(pDiffr)) - offset;
        end
        targetData = CreateFrequencyNBands(validationData, fvec, nBands);
        targetData = max(-128, min(128, targetData));
        % trainingData = max(-128, min(128, trainingData)) / 128;
        if (saveValidation)
            save([cd filesep savePath '.mat'], 'trainingData', 'targetData', 'validationData', 'geometry');
        else
            save([cd filesep savePath '.mat'], 'trainingData', 'targetData');
        end
        % save([cd filesep savePath2 '.mat'], 'fc', 'gain', 'blendExpn', 'Q', 'udfaTerms');
        toc
    end

    % Transform Fc
    trainingData(1:numNNInputs / 2,:) = log10(trainingData(1:numNNInputs / 2,:));

    % Transform gain
    gainSign = sign(trainingData(numNNInputs / 2:numNNInputs,:));
    trainingData(numNNInputs / 2 + 1:numNNInputs,:) = log10(1+abs(trainingData(numNNInputs / 2 + 1:numNNInputs,:)));

    % Clip max value
    trainingData = min(trainingData, 10);
    trainingData(numNNInputs / 2:numNNInputs,:) = gainSign .* trainingData(numNNInputs / 2:numNNInputs,:);
end