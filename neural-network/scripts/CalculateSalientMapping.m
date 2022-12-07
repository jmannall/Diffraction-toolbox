%% Biquad loss function for training neural networks

function gradients = CalculateSalientMapping(net, inputData, targetData, lossFunc)

    output = predict(net, inputData);
    [numOutputs, numInputs] = size(output);

    for i = 1:numOutputs
        gradients{i} = dlgradient(sum(output(i,:)), inputData);
    end
end