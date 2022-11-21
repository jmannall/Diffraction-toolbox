function [net, losses] = ParTrainNeuralNetwork(net, trainingData, targetData, epoch, miniBatchSize, numIterationsPerEpoch, lossFunc, x)
    
    learnRate = x.lR;
    gradDecay = x.gD;
    sqGradDecay = x.sGD;
    maxGradient = x.mG;

    averageGrad = [];
    averageSqGrad = [];
    losses = zeros(1, numIterationsPerEpoch);

    % Shuffle data.
    idx = randperm(size(trainingData, 2));
    targetData = targetData(:,idx);
    trainingData = trainingData(:,idx);
    for i = 1:numIterationsPerEpoch
        iteration = (epoch - 1) * numIterationsPerEpoch + i;

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
        
        loss = double(loss);
        losses(i) = loss;
    end
end