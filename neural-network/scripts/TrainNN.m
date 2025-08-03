function [net, losses] = TrainNN(net, hP, tP, nP, losses)

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
    thisIterationLosses = zeros(1, numIterationsPerEpoch);

    iteration = 0;
    start = tic;
    i = 1;

    lines = CreateAnimatedLinePlot();

    %numLayers = length(net.Learnables.Value);
    %numNodes = numel(net.Learnables.Value{2});

    % File save data
    fileName = CreateNNIdx(hP, tP, nP);
    tempDir = 'tempNN';
    CheckFileDir(tempDir);
    savePath = [tempDir filesep fileName];
    loadPath = [savePath '.mat'];
    backupRate = 5;
 
    % Load testing data
    [testInputData, testTargetData] = tP.dataFunc('TestData');
    %[testInputData, testTargetData, DC] = tP.dataFunc('TestData');
    %testTargetData = [DC; testTargetData];
    testInputData = dlarray(single(testInputData), "CB");

    %reductionPoints = [100 300 450];
    %reductionPoints = [40 75 90];
    reductionPoints = [15 40 75 90];

    % Check if restarting training
    restart = exist([cd filesep loadPath], "file");
    if restart == 2
        load(loadPath, "net", "losses", "hP", "tP", "nP", "iteration", "i")
        disp('Restart training')
    elseif losses.iteration(1) > 0
        i = find(losses.epoch == 1000, 1);
        iteration = find(losses.iteration == 0, 1);
        disp('Continue training')
    else
        disp('Start training')
    end

    tic
    for epoch = i:tP.numEpochs
        [trainingData, targetData] = tP.dataFunc(epoch);
        %[trainingData, targetData, DCTargetData] = tP.dataFunc(epoch);
        targetData = dlarray(targetData);
        % Shuffle data.
        idx = randperm(size(trainingData, 2));
        targetData = targetData(:,idx);
        trainingData = trainingData(:,idx);
        %DCTargetData = DCTargetData(idx);
        %targetData = [DCTargetData; targetData];
        for j = 1:numIterationsPerEpoch
            iteration = iteration + 1;
    
            % Create input and target variables
            [X, T] = ProcessNNBatchInputData(trainingData, targetData, tP.miniBatchSize, j);
    
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
            [net,averageGrad,averageSqGrad] = adamupdate(net, gradients, averageGrad, averageSqGrad, iteration, hP.learnRate, tP.gradDecay, tP.sqGradDecay);
            
            % Store losses
            idx = numIterationsPerEpoch * (epoch - 1) + j;
            loss = double(loss);
            losses.iteration(idx) = loss;
            thisIterationLosses(j) = loss;
        end
        losses.test(epoch) = tP.testFunc(net, testInputData, testTargetData);
        losses.epoch(epoch) = mean(thisIterationLosses);
        if sum(ismember(epoch, reductionPoints)) > 0
            hP.learnRate = hP.learnRate / 10;
        end

        UpdateNNAnimatedLinePlot(lines, losses, numIterationsPerEpoch, epoch, start)

        % Backup result every 5 epochs
        if mod(epoch, backupRate) == 0
            i = epoch + 1;
            worker = getCurrentWorker;
            if ~isempty(worker)
                %disp(['Save backup: ', num2str(worker.ProcessId)])
            else
                %disp('Save backup')
            end
            save(savePath, "net", "losses", "hP", "tP", "nP", "iteration", "i", '-v7.3')
        end
        if epoch > 9 && losses.test(epoch) > 1e3
            disp(['Early stop criteria reached. Epoch: ' num2str(epoch)])
            nP.savePath = [nP.savePath '_INCOMPLETE'];
            return
        end
        if isnan(loss)
            disp(['End training early. Epoch: ' num2str(epoch)])
            nP.savePath = [nP.savePath '_INCOMPLETE'];
            return
        end
    end
    toc

    disp('Save trained network')
    save(nP.savePath, "net", "losses", "hP", "tP", "nP", '-v7.3')

    % Delete temp backup
    delete(loadPath)
end