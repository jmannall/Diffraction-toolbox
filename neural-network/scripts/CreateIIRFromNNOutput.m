%% Create IIR filter parameters from a neural network output

function [z, p, k] = CreateIIRFromNNOutput(output, filterOrder)
    z = output(1:filterOrder,:);
    p = output(filterOrder + 1:2 * filterOrder,:);

    epsilon = 1e-8;
    magZ = abs(z);
    magP = abs(p);
    z = (1 - epsilon) .* z .* tanh(magZ) ./ (magZ + epsilon);
    p = (1 - epsilon) .* p .* tanh(magP) ./ (magP + epsilon);
    k = Sigmoid(output(end,:));
end