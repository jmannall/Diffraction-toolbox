function [b, a] = EDlpf(fc, T)
    omega = 2 * pi * fc;
    K = omega * T;
    
    b = [K; K];
    a = [K + 2; K - 2];
end