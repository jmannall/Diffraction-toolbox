function [b, a] = LowPassCoefficients(fc, fs, k)

    T = 1 / fs;
    omega = 2 * pi * fc;
    K = omega * T;
    
    b = [K; K];
    a = [K + 2; K - 2];
    if nargin > 2
        b = b .* k;
    end
    [x, z] = size(b);
    b = reshape(b, x, 1, z);
    a = reshape(a, x, 1, z);
end