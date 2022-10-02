%% Leaky ReLU activation function

function output = Relu(x)
    output = max(0.2 * x, x);
end