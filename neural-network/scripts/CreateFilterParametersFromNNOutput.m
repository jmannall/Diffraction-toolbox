function [lpFc, hsFc, G, k] = CreateFilterParametersFromNNOutput(output)

    lpFc = 10e4 * Sigmoid(output(1,:));
    hsFc = 14e4 * Sigmoid(output(2,:));
    G = output(3,:);
    k = Sigmoid(output(4,:));
end