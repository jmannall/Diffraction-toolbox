%% Derivative of sigmoid activation function

function output = Delsigmoid(x)
    output = Sigmoid(x) .* (1 - Sigmoid(x));
end