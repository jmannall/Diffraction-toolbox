function [net, epochLosses, losses] = TrainNeuralNetwork(net, trainingData, targetData, numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x, dataFunc)
    
    learnRate = x.lR;
    gradDecay = x.gD;
    sqGradDecay = x.sGD;
    maxGradient = x.mG;

    averageGrad = [];
    averageSqGrad = [];
    losses = zeros(1, numEpochs * numIterationsPerEpoch);
    iterationLosses = zeros(1, numIterationsPerEpoch);
    epochLosses = zeros(1, numEpochs);
    iteration = 0;
    start = tic;
    count = 1;
    i = 1;

    [lineIterationLoss, lineEpochLoss] = CreateAnimatedLinePlot();
    oldNet = [];

    if nargin < 9
        idx = DataHash({trainingData, extractdata(targetData), numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x, dataFunc});
    else
        idx = DataHash({trainingData, extractdata(targetData), numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x});
    end

    filePath = 'tempNN';
    CheckFileDir(filePath);
    fileName = num2str(idx);
    savePath = [filePath filesep fileName];
    loadPath = [savePath '.mat'];

    backupRate = 5;
    restart = exist([cd filesep loadPath], "file");

    if restart == 2
        load(loadPath, "net", "losses", "epochLosses", "iteration", "i")
        disp('Restart training')
    else
        disp('Start training')
    end

    tic
    for epoch = i:numEpochs
        if nargin > 8 && epoch > 1
            [trainingData, targetData] = dataFunc(epoch);
        end
        % Shuffle data.
        idx = randperm(size(trainingData, 2));
        targetData = targetData(:,idx);
        trainingData = trainingData(:,idx);
        for i = 1:numIterationsPerEpoch
            iteration = iteration + 1;
    
            [X, T] = ProcessNNBatchInputData(trainingData, targetData, miniBatchSize, i);
    
            % Evaluate the model loss and gradients using dlfeval and the
            % modelLoss function.input
            [loss, state, gradients] = dlfeval(lossFunc,net,X,T);
            % Update normalisation paramters
%             x = state.Value(2);
%             y = extractdata(x{1})
            isVariance = strcmp(state.Parameter, "TrainedVariance");
            state.Value(isVariance) = cellfun(@(x) max(x, 1e-10), state.Value(isVariance), 'UniformOutput', false);
            net.State = state;
    
            % Update the network parameters using the Adam optimizer.
            % Clip gradients
            gradients = ClipGradients(gradients, maxGradient);
            [net,averageGrad,averageSqGrad] = adamupdate(net,gradients,averageGrad,averageSqGrad,iteration,learnRate,gradDecay,sqGradDecay);
            
            idx = numIterationsPerEpoch * (epoch - 1) + i;
            loss = double(loss);
            losses(idx) = loss;
            iterationLosses(i) = loss;
        end
        if isnan(loss)
            % Revert results to last epoch and end training
            lastEpoch = max(epoch - 1, 1);
            epochLosses = epochLosses(1:lastEpoch);
            losses = losses(1:lastEpoch * numIterationsPerEpoch);
            net = oldNet;
            if epoch > backupRate
                delete(loadPath)
            end
            break
        end
        epochLosses(epoch) = mean(iterationLosses);
        if count > 10
            idx = max(epoch - 10, 1):epoch - 1;
            runningMean = mean(epochLosses(idx));
            improvement = runningMean - epochLosses(epoch);
            count = 1;
        else
            improvement = 1;
            count = count + 1;
        end
        if improvement < 0
            learnRate = learnRate / 2;
        end

        UpdateNNAnimatedLinePlot(lineIterationLoss, lineEpochLoss, losses, numIterationsPerEpoch, epoch, start)
        oldNet = net;

        % Backup result every 5 epochs
        if mod(epoch, backupRate) == 0
            i = epoch;
            save(savePath, "net", "losses", "epochLosses", "iteration", "i", '-v7.3')
        end
    end
    delete(loadPath)
    toc
end