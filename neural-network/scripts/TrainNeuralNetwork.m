function [net, epochLosses, losses] = TrainNeuralNetwork(net, trainingData, targetData, numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x, restartFunc, dataFunc)
    
    learnRate = x.lR;
    gradDecay = x.gD;
    sqGradDecay = x.sGD;
    maxGradient = x.mG;

    worker = getCurrentWorker;
    if ~isempty(worker)
        disp(['Worker: ', num2str(worker.ProcessId), ', Network: ' num2str(x.nL), '-', num2str(x.hL)]);
    else
        disp(['Network: ' num2str(x.nL), '-', num2str(x.hL)]);
    end
%     gradDecay = 0.9;
%     sqGradDecay = 0.99;
%     maxGradient = 1.5;

    averageGrad = [];
    averageSqGrad = [];
    losses = zeros(1, numEpochs * numIterationsPerEpoch);
    iterationLosses = zeros(1, numIterationsPerEpoch);
    epochLosses = zeros(1, numEpochs);
    iteration = 0;
    start = tic;
    i = 1;
    count = 0;
    firstEpoch = true;

    [lineIterationLoss, lineEpochLoss] = CreateAnimatedLinePlot();
    oldNet = [];

    numLayers = length(net.Learnables.Value);
    numNodes = numel(net.Learnables.Value{2});
    if nargin < 9
        idx = DataHash({numLayers, numNodes, trainingData, targetData, numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x, restartFunc, dataFunc});
    else
        idx = DataHash({numLayers, numNodes, trainingData, targetData, numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x, restartFunc});
    end

    filePath = 'tempNN';
    CheckFileDir(filePath);
    fileName = num2str(idx);
    savePath = [filePath filesep fileName];
    loadPath = [savePath '.mat'];

    backupRate = 10;
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
        targetData = dlarray(targetData);
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
        epochLosses(epoch) = mean(iterationLosses);
        improvement = epochLosses(max(1, epoch - 1)) - epochLosses(epoch);
        if epoch > 1 && (isnan(loss) || improvement < -100)
            % Revert results to last epoch and end training
            if firstEpoch
                if ~isempty(worker)
                    disp(['Reinitialise net: ', num2str(epoch), '_', num2str(worker.ProcessId)])
                else
                    disp(['Reinitialise net: ', num2str(epoch)])
                end
                net = restartFunc();
            else
                firstEpoch = false;
                if ~isempty(worker)
                    disp(['Return to previous epoch net: ', num2str(epoch), '_', num2str(worker.ProcessId)])
                else
                    disp(['Return to previous epoch net: ', num2str(epoch)])
                end
                net = oldNet;
            end
            count = count + 1;
        else
            firstEpoch = false;
            count = 0;
        end
        if count == 5
            disp(['Reduce learn rate early: ', num2str(epoch)])
            learnRate = learnRate / 10;
        end
        if count > 10
            lastEpoch = max(epoch - 1, 1);
            epochLosses = epochLosses(1:lastEpoch);
            losses = losses(1:lastEpoch * numIterationsPerEpoch);
            disp(['End training early: ', num2str(epoch)])
            break
        end
        if epoch == 200 || epoch == 400
            count = 0;
            learnRate = learnRate / 10;
        end

        UpdateNNAnimatedLinePlot(lineIterationLoss, lineEpochLoss, losses, numIterationsPerEpoch, epoch, start)
        oldNet = net;

        % Backup result every 10 epochs
        if mod(epoch, backupRate) == 0
            i = epoch;
            worker = getCurrentWorker;
            if ~isempty(worker)
                disp(['Save backup: ', num2str(worker.ProcessId)])
            else
                disp('Save backup')
            end
            save(savePath, "net", "losses", "epochLosses", "iteration", "i", '-v7.3')
        end
    end
    if epoch > backupRate
        delete(loadPath)
    end    
    toc
end