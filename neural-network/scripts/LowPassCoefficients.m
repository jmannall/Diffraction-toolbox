function [b, a] = LowPassCoefficients(fc, fs, k)

    T = 1 / fs;
    omega = 2 * pi * fc;
    K = omega * T;
    
    b(:,1,:) = [K; K];
    a(:,1,:) = [K + 2; K - 2];
    if nargin > 2
        b(:,1,:) = squeeze(b(:,1,:)) .* k;
    end
end