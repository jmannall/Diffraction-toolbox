function BayesoptNeuralNetwork(lossFunc, networkSize, numOutputs, controlparameters)
    
    maxLayers = networkSize / 2 ^ 2;
    
    learnRate = optimizableVariable('lR',[1e-5 1e-1],'Type','real','Transform','log');
    gradDecay = optimizableVariable('gD',[0.5 1],'Type','real');
    sqGradDecay = optimizableVariable('sGD',[0.8 1],'Type','real');
    maxGradient = optimizableVariable('mG',[0.5 10],'Type','real','Transform','log');
    numLayers = optimizableVariable('nL',[2 min(20, maxLayers)],'Type','integer');
    
    epochSize = 20e3;

    saveDir = 'runningFiles';
    CheckFileDir(saveDir);
    disp(['Save path ', saveDir]);
    saveSeed = [saveDir, '\Seed.mat'];
    saveResult = [saveDir, '\BayesoptResults.mat'];
    restarting = isfile(saveSeed);
    if restarting
        load(saveSeed, 'seed')
    else
        seed = round(1e9 * rand(1));
        save(saveSeed, 'seed');
    end
    rng(seed)
    
    disp('Create Taining Data')
    [~, ~, ~, ~, ~, idx, saveData] = CreateBtmTrainingData(epochSize, controlparameters);
    dataFunc = @() CreateBtmTrainingData(epochSize, controlparameters, idx);

    numEpochs = 100;
    func = @(x)CreateBTMNeuralNetwork(x, lossFunc, dataFunc, networkSize, numOutputs, numEpochs);
    restarting = isfile(saveResult);
    if restarting
        disp('Restart Optimisation')
        load(saveResult, 'BayesoptResults');
        result = resume(BayesoptResults, 'SaveFileName', saveResult, 'OutputFcn', @saveToFile );
    else
        disp('Start Optimisation')
        result = bayesopt(func, [learnRate, gradDecay, sqGradDecay, maxGradient, numLayers], ...
        'UseParallel', true, 'MaxTime', 86400, 'MaxObjectiveEvaluations', 30, ...
        'SaveFileName', saveResult, 'OutputFcn', @saveToFile);
    end
    
    xObs = result.XAtMinObjective;
    xEst = result.XAtMinEstimatedObjective;
    
    %% Test
    
    [lossObs, netObs] = func(xObs);
    [lossEst, netEst] = func(xEst);
    
    %% Save
    
    save(['BayesoptResult_NNSize_', num2str(networkSize), '.mat'], "result", "xObs", "xEst", "lossObs", "netObs", "lossEst", "netEst")
    
    %% Delete data save
    delete([saveData, '.mat'])
    delete(saveSeed)
    delete(saveResult)
end