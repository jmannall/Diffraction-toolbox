%% Biquad loss function for training neural networks

function [loss, state, gradients] = NNBiquadLoss(net, inputData, targetData, numBiquads, numFreq, fs, fidx, training)

    if training
        [output, state] = forward(net, inputData);
    else
        output = predict(net, inputData);
        state = 0;
    end
    
    loss = BiquadLoss(output, targetData, numBiquads, numFreq, fs, fidx);
    gradients = dlgradient(loss, net.Learnables);
end