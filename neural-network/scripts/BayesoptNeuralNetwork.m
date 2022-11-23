function BayesoptNeuralNetwork(lossFunc, networkSize, numOutputs, controlparameters)
    
    maxLayers = networkSize / 2 ^ 2;
    
    learnRate = optimizableVariable('lR',[1e-5 1e-1],'Type','real','Transform','log');
    gradDecay = optimizableVariable('gD',[0.5 1],'Type','real');
    sqGradDecay = optimizableVariable('sGD',[0.8 1],'Type','real');
    maxGradient = optimizableVariable('mG',[0.5 10],'Type','real','Transform','log');
    numLayers = optimizableVariable('nL',[2 min(20, maxLayers)],'Type','integer');
    
    epochSize = 200;

    saveDir = 'runningFiles';
    CheckFileDir(saveDir);
    saveResult = [saveDir, '\BayesoptResults.mat'];
    
    disp('Create Taining Data')
    [~, ~, ~, ~, ~, idx, saveData] = CreateBtmTrainingData(epochSize, controlparameters);
    dataFunc = @() CreateBtmTrainingData(epochSize, controlparameters, idx);

    numEpochs = 5;
    func = @(x)CreateBTMNeuralNetwork(x, lossFunc, dataFunc, networkSize, numOutputs, numEpochs);
    restarting = isfile(saveResult);
    if restarting
        disp('Restart Optimisation')
        load(saveResult, 'BayesoptResults');
        result = resume(BayesoptResults, 'SaveFileName', saveResult, 'OutputFcn', @saveToFile );
    else
        disp('Start Optimisation')
        result = bayesopt(func, [learnRate, gradDecay, sqGradDecay, maxGradient, numLayers], ...
        'UseParallel', true, 'MaxTime', 86400, 'MaxObjectiveEvaluations', 20, ...
        'SaveFileName', saveResult, 'OutputFcn', @saveToFile);
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
    
    save(['BayesoptResult_Size_', num2str(networkSize), '_Filter_', controlparameters.filterType, '.mat'], "result", "xObs", "xEst", "lossObs", "netObs", "lossEst", "netEst")
    
    %% Delete data save
    delete([saveData, '.mat'])
    delete(saveResult)
end