function [a, b] = NNOutputToFilterCoefficients(out)
    [z, p, k] = TransformNNOutput(out);
    b = k * [1, -sum(z), prod(z)];
    a = [1, -sum(p), prod(p)];
end