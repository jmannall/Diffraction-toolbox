function [b, a] = EDhsh(fc, G, T)
    Bpi = 10 ^ (G / 20);
    omega = 2 * pi * fc;
    K = omega * T;
    
    b = [K + 2 * Bpi; K - 2 * Bpi];
    a = [K + 2; K - 2];
end