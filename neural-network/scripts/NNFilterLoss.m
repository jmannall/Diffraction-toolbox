%% Biquad loss function for training neural networks

function [loss, state, gradients, predictionN, prediction] = NNFilterLoss(net, inputData, targetData, lossFunc, training)

    if training
        [output, state] = forward(net, inputData);
    else
        output = predict(net, inputData);
        state = 0;
    end
    
    [loss, predictionN, prediction] = lossFunc(output, targetData);

    if training
        gradients = dlgradient(loss, net.Learnables);
    else
        loss = extractdata(loss);
        predictionN = extractdata(predictionN);
        prediction = extractdata(prediction);
        gradients = 0;
    end
end