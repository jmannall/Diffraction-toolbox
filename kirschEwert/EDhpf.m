function [b, a] = EDhpf(fc, T)
    omega = 2 * pi * fc;
    K = omega * T;
    
    b = [2; -2];
    a = [K + 2; K - 2];
end