%% Create IIR filter parameters from a neural network output

function [z, p, k] = CreateIIRFromNNOutput(output, numIIRFilters)
    z = output(1:numIIRFilters,:);
    p = output(numIIRFilters + 1:2 * numIIRFilters,:);

    epsilon = 1e-10;
    z = (1 - epsilon) .* z .* tanh(z) ./ (z + epsilon);
    p = (1 - epsilon) .* p .* tanh(p) ./ (p + epsilon);
    k = Sigmoid(output(end,:));
end