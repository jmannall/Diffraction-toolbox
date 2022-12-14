function [b, a] = HighShelfCoefficients(fc, G, fs, k)

    T = 1 / fs;
    omega = 2 * pi * fc;
    K = omega * T;
    B0 = 10 .^ (G / 20);

    b = [K + 2 * B0; K - 2 * B0];
    a = [K + 2; K - 2];
    if nargin > 3
        b = b .* k;
    end
    [x, z] = size(b);
    b = reshape(b, x, 1, z);
    a = reshape(a, x, 1, z);
end
