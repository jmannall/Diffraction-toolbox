function [net, epochLosses, losses] = TrainNeuralNetwork(net, trainingData, targetData, numEpochs, miniBatchSize, numIterationsPerEpoch, lossFunc, x, dataFunc)
    
    learnRate = x.lR;
    gradDecay = x.gD;
    sqGradDecay = x.sGD;
    maxGradient = x.mG;

    averageGrad = [];
    averageSqGrad = [];
    losses = zeros(1, numEpochs * numIterationsPerEpoch);
    iterationLosses = zeros(1, numIterationsPerEpoch);
    epochLosses = zeros(1, numIterationsPerEpoch);
    iteration = 0;
    start = tic;

    [lineIterationLoss, lineEpochLoss] = CreateAnimatedLinePlot();

    disp('Start training')
    tic
    for epoch = 1:numEpochs
        if nargin > 8 && epoch > 1
            [trainingData, targetData] = dataFunc();
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
        epochLosses(epoch) = mean(iterationLosses);
        UpdateNNAnimatedLinePlot(lineIterationLoss, lineEpochLoss, losses, numIterationsPerEpoch, epoch, start)
    end
    toc
end