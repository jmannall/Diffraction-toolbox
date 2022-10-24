function [loss, state, gradients] = MeanSquaredError(net, inputData, targetData)

    [output, state] = forward(net, inputData);
    
    loss = sum((output - targetData).^2 / (numel(targetData)), 'all');
    gradients = dlgradient(loss, net.Learnables);
end