function [z, p, k] = TransformNNOutput(output)
    z = output(1:2,:);
    p = output(3:4,:);

    epsilon = 1e-8;
    magZ = abs(z);
    magP = abs(p);
    z = (1 - epsilon) .* z .* tanh(magZ) ./ (magZ + epsilon);
    p = (1 - epsilon) .* p .* tanh(magP) ./ (magP + epsilon);
    k = Sigmoid(output(5,:));
end