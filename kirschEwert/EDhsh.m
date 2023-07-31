function [b, a] = EDhsh(fc, G, T)
    omega = 2 * pi * fc;
    omegaD = cot(omega * T / 2);
    
    b = [1 + sqrt(G) * omegaD; 1 - sqrt(G) * omegaD];
    a = [1 + omegaD / sqrt(G); 1 - omegaD / sqrt(G)];

    b = b / a(1);
    a = a / a(1);
end