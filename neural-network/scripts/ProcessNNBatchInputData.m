function [X, T] = ProcessNNBatchInputData(trainingData, targetData, miniBatchSize, i) 

    % Read mini-batch of data and convert the labels to dummy
    % variables.
    idx = (i-1)*miniBatchSize+1:i*miniBatchSize;
    X = trainingData(:,idx);
    T = targetData(:,idx);

    % Convert mini-batch of data to a dlarray.
    X = dlarray(single(X), "CB");

    % If training on a GPU, then convert data to a gpuArray.
    if canUseGPU
        X = gpuArray(X);
    end
end