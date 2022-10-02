%% Derivative of leaky ReLU activation function

function output = Delrelu(x)
    output = ones(size(x));
    idx = x < 0;
    output(idx) = 0.2;
end