%% MSE error loss function

function [loss, delloss] = MSELoss(input, target)
    loss = 0.5 * (input - target).^2;
    delloss = dlgradient(sum(loss), input);
end