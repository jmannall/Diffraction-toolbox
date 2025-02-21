function TrainUDFA_NN(numLayers, size, learnRate)

    close all

    if learnRate > 1
        learnRate = 1 / learnRate;
    end

    % Control parameters
    controlparameters.fs = 48e3;
    controlparameters.nfft = 8192;
    controlparameters.filterOrder = 2;
    controlparameters.nBands = 8;
    controlparameters.fvec = controlparameters.fs/controlparameters.nfft*[0:controlparameters.nfft/2-1];
    controlparameters.fidx = CreateFidx(controlparameters.fvec, controlparameters.nBands);
    controlparameters.numNNInputs = 8;

    filterFunc = @(output, target) IIRFilterLoss(output, target, controlparameters);
    
    % Training parameters
    numEpochs = 50;
    epochSize = 5e3; %20e3
    miniBatchSize = 128;

    gradDecay = 0.9;
    sqGradDecay = 0.99;
    maxGrad = 0.9;

    lossFunc = @(net, trainingData, targetData) NNFilterLoss(net, trainingData, targetData, filterFunc, true);
    dataFunc = @(idx) CreateUDFA_NNTrainingData(epochSize, controlparameters, false, idx);
    testFunc = @(net, trainingData, targetData) NNFilterLoss(net, trainingData, targetData, filterFunc, false);

    % Network parameters
    numInputs = controlparameters.numNNInputs;
    numOutputs = 2 * controlparameters.filterOrder + 1;
    alpha = 0.2;
    
    % NN complexity
    gx = 5;
    CcPZTransform = 6;
    numPZ = numOutputs - 1;
    CcKTransform = 1;
    CcOutputTransform = numPZ * CcPZTransform + CcKTransform;
    CcIIR = 6;
    
    a = numLayers - 1;
    b = numLayers + numLayers * gx + numInputs + numOutputs;
    c = CcIIR + CcOutputTransform + numOutputs - size;
    
    hiddenLayerSize = round((-b + sqrt(b .^ 2 - 4 .* c * a)) ./ (2 * a));
    
    % Save paths
    rootDir = 'NNSaves';
    CheckFileDir(rootDir)
    CheckFileDir([rootDir filesep saveDir])
    idx = [num2str(numLayers), '_', num2str(hiddenLayerSize), '_', num2str(learnRate)];
    idx = erase(idx, '.');
    saveFile = ['NN-', idx];
    savePath = [rootDir filesep saveDir filesep saveFile];

    rng shuffle
    seed = rng().Seed;

    % Parameter structs
    hyperParameters = struct('learnRate', learnRate, 'numLayers', numLayers, 'hiddenLayerSize', hiddenLayerSize);
    trainingParameters = struct('lossFunc', lossFunc, 'dataFunc', dataFunc, 'testFunc', testFunc, 'numEpochs', numEpochs, 'epochSize', epochSize, ...
        'miniBatchSize', miniBatchSize, 'gradDecay', gradDecay, 'sqGradDecay', sqGradDecay, 'maxGrad', maxGrad);
    networkParamters = struct('numInputs', numInputs, 'numOutputs', numOutputs, 'alpha', alpha, 'saveDir', saveDir, 'savePath', savePath, 'seed', seed);

    complete1 = exist([cd filesep savePath '.mat'], "file");
    complete2 = exist([cd filesep savePath '_INCOMPLETE.mat'], "file");

    if complete2 == 2
        disp('Already Trained')
    elseif complete1 == 2
        load([cd filesep savePath '.mat'], 'net', 'hP', 'tP', 'nP', 'losses')
        if tP.numEpochs < numEpochs
            tP = trainingParameters;
            [loss, net] = CreateNN(hP, tP, nP, net, losses);
        else
            disp('Already Trained')
        end
    else
        % Train the network
        [loss, net] = CreateNN(hyperParameters, trainingParameters, networkParamters);
    end
end
