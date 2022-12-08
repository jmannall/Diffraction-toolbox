function BayesoptNeuralNetwork(lossFunc, networkSizeInput, numOutputs, controlparameters)
    
    maxLayers = networkSizeInput / 2 ^ 2;
    
    learnRate = optimizableVariable('lR',[1e-5 1e-2],'Type','real','Transform','log');
    gradDecay = optimizableVariable('gD',[0.5 1],'Type','real');
    sqGradDecay = optimizableVariable('sGD',[0.8 1],'Type','real');
    maxGradient = optimizableVariable('mG',[0.5 10],'Type','real','Transform','log');
    numLayers = optimizableVariable('nL',[2 min(20, maxLayers)],'Type','integer');
    networkSize = optimizableVariable('nS',[5000 networkSizeInput],'Type','real','Transform','log');
    
    epochSize = 20e3;

    saveDir = 'runningFiles';
    CheckFileDir(saveDir);
    saveResult = [saveDir, filesep, 'BayesoptResults.mat'];
    
    disp('Create Training Data')
    % [~, ~, ~, ~, ~, idx, saveData] = CreateBtmTrainingData(epochSize, controlparameters, epochSize);
    dataFunc = @(idx) CreateBtmTrainingData(epochSize, controlparameters, idx);

    numEpochs = 100;
    func = @(x)CreateBTMNeuralNetwork(x, lossFunc, dataFunc, networkSizeInput, numOutputs, numEpochs, controlparameters.filterType, true);
    restarting = isfile(saveResult);
    if restarting
        disp('Restart Optimisation')
        load(saveResult, 'BayesoptResults');
        result = resume(BayesoptResults, 'SaveFileName', saveResult, 'OutputFcn', @saveToFile );
    else
        disp('Start Optimisation')
        result = bayesopt(func, [learnRate, gradDecay, sqGradDecay, maxGradient, numLayers, networkSize], ...
        'UseParallel', true, 'MaxTime', 86400, 'MaxObjectiveEvaluations', 100, ...
        'SaveFileName', saveResult, 'OutputFcn', @saveToFile, 'AcquisitionFunctionName', 'expected-improvement-plus');
    end
    
    xObs = result.XAtMinObjective;
    xEst = result.XAtMinEstimatedObjective;
    
    %% Test
    
    x = [xObs; xEst];
    parfor i = 1:2
        [loss(i), net(i)] = func(x(i,:));
    end

    lossObs = loss(1);
    lossEst = loss(2);
    netObs = net(1);
    netEst = net(2);

    %% Save
    
    saveDir = 'bayesoptResults';
    CheckFileDir(saveDir)
    save([saveDir, filesep, 'BayesoptResult_Size_', num2str(networkSize), '_Filter_', controlparameters.filterType, '.mat'], "result", "xObs", "xEst", "lossObs", "netObs", "lossEst", "netEst")
    
    %% Delete data save
    % delete([saveData, '.mat'])
    delete(saveResult)
end