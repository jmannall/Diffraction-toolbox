%% Biquad loss function for training neural networks

function [loss, state, gradients] = NNBiquadLoss(net, inputData, targetData, numBiquads, numFreq, fidx, training)

    if training
        [output, state] = forward(net, inputData);
    else
        output = predict(net, inputData);
        state = 0;
    end
    
    [loss, gradients] = BiquadLoss(output, targetData, numBiquads, numFreq, fidx);
end