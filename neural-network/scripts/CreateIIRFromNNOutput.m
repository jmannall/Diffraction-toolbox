%% Create IIR filter parameters from a neural network output

function [z, p, k] = CreateIIRFromNNOutput(output, numIIRFilters)
    z = output(1:numIIRFilters,:);
    p = output(numIIRFilters + 1:2 * numIIRFilters,:);

    epsilon = 1e-8;
    magZ = abs(z);
    magP = abs(p);
    z = (1 - epsilon) .* z .* tanh(magZ) ./ (magZ + epsilon);
    p = (1 - epsilon) .* p .* tanh(magP) ./ (magP + epsilon);
    k = 100 * Sigmoid(output(end,:));
end