%% Biquad loss function for training neural networks

function [loss, state, gradients] = NNFilterLoss(net, inputData, targetData, lossFunc, training)

    if training
        [output, state] = forward(net, inputData);
    else
        output = predict(net, inputData);
        state = 0;
    end
    
    loss = lossFunc(output, targetData);
    gradients = dlgradient(loss, net.Learnables);
end