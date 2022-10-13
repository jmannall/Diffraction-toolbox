%% Create biquad filter parameters from a neural network output

function [zR, zI, pR, pI, k] = CreateZPKFromNNOutput(output, numBiquads)
    [zR, zI] = MinPhaseTransform(output(1:numBiquads,:), output(numBiquads + 1:2 * numBiquads,:));
    [pR, pI] = MinPhaseTransform(output(2 * numBiquads + 1:3 * numBiquads,:), output(3 * numBiquads + 1:4 * numBiquads,:));
    k = Sigmoid(output(end,:));
end