function [net, losses] = TrainNN(net, hP, tP, nP)

    worker = getCurrentWorker;
    name = ['Network: ' num2str(hP.numLayers), '-', num2str(hP.hiddenLayerSize), '-' num2str(hP.learnRate)];
    if ~isempty(worker)
        disp(['Worker: ', num2str(worker.ProcessId), ', ', name]);
    else
        disp(name);
    end

    % Initialise gradient variables
    averageGrad = [];
    averageSqGrad = [];

    % Initialise losss variables
    numIterationsPerEpoch = floor(tP.epochSize./tP.miniBatchSize);
    losses.iteration = zeros(1, tP.numEpochs * numIterationsPerEpoch);
    [losses.epoch, losses.test] = deal(zeros(1, tP.numEpochs));
    thisIterationLosses = zeros(1, numIterationsPerEpoch);

    iteration = 0;
    start = tic;
    i = 1;

    [lineIterationLoss, lineEpochLoss] = CreateAnimatedLinePlot();

    %numLayers = length(net.Learnables.Value);
    %numNodes = numel(net.Learnables.Value{2});

    % File save data
    idx = DataHash({hP, tP, nP});

    tempDir = 'tempNN';
    CheckFileDir(tempDir);
    fileName = num2str(idx);
    savePath = [tempDir filesep fileName];
    loadPath = [savePath '.mat'];
    backupRate = 10;

    % Load testing data
    [testInputData, testTargetData] = tP.dataFunc('TestData');
    testInputData = dlarray(single(testInputData), "CB");
    learnRate = hP.learnRate;

    % Check if restarting training
    restart = exist([cd filesep loadPath], "file");
    if restart == 2
        load(loadPath, "net", "iterationLosses", "epochLosses", "iteration", "i")
        disp('Restart training')
    else
        disp('Start training')
    end

    tic
    for epoch = i:tP.numEpochs
        [trainingData, targetData] = tP.dataFunc(epoch);
        targetData = dlarray(targetData);
        % Shuffle data.
        idx = randperm(size(trainingData, 2));
        targetData = targetData(:,idx);
        trainingData = trainingData(:,idx);
        for i = 1:numIterationsPerEpoch
            iteration = iteration + 1;
    
            % Create input and target variables
            [X, T] = ProcessNNBatchInputData(trainingData, targetData, tP.miniBatchSize, i);
    
            % Evaluate the model loss and gradients using dlfeval and the
            % modelLoss function.input
            [loss, state, gradients] = dlfeval(tP.lossFunc, net, X, T);

            % Update normalisation paramters
            isVariance = strcmp(state.Parameter, "TrainedVariance");
            state.Value(isVariance) = cellfun(@(x) max(x, 1e-10), state.Value(isVariance), 'UniformOutput', false);
            net.State = state;
    
            % Clip gradients
            gradients = ClipGradients(gradients, tP.maxGrad);
            % Update the network parameters using the Adam optimizer.
            [net,averageGrad,averageSqGrad] = adamupdate(net, gradients, averageGrad, averageSqGrad, iteration, learnRate, tP.gradDecay, tP.sqGradDecay);
            
            % Store losses
            idx = numIterationsPerEpoch * (epoch - 1) + i;
            loss = double(loss);
            losses.iteration(idx) = loss;
            thisIterationLosses(i) = loss;
        end
        losses.test(epoch) = tP.testFunc(net, testInputData, testTargetData);
        losses.epoch(epoch) = mean(thisIterationLosses);
        if epoch == 200 || epoch == 400
            learnRate = learnRate / 10;
        end

        UpdateNNAnimatedLinePlot(lineIterationLoss, lineEpochLoss, losses, numIterationsPerEpoch, epoch, start)

        % Backup result every 10 epochs
        if mod(epoch, backupRate) == 0
            i = epoch;
            worker = getCurrentWorker;
            if ~isempty(worker)
                disp(['Save backup: ', num2str(worker.ProcessId)])
            else
                disp('Save backup')
            end
            save(savePath, "net", "iterationLosses", "epochLosses", "iteration", "i", '-v7.3')
        end
        if isnan(loss)
            disp(['End training early. Epoch: ' num2str(epoch)])
            nP.savePath = [nP.savePath '_INCOMPLETE'];
            break
        end
    end
    toc

    disp('Save trained network')
    save(nP.savePath, "net", "losses", "hP", "tP", "nP", '-v7.3')

    % Delete temp backup
    delete(loadPath)
end