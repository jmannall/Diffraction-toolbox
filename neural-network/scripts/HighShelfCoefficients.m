function [b, a] = HighShelfCoefficients(fc, G, fs, k)

    T = 1 / fs;
    omega = 2 * pi * fc;
    K = omega * T;
    B0 = 10 ^ (G / 20);

    b(:,1,:) = [K + 2 * B0; K - 2 * B0];
    a(:,1,:) = [K + 2; K - 2];
    if nargin > 3
        b(:,1,:) = squeeze(b(:,1,:)) .* k;
    end
end
