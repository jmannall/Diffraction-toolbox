function [b, a] = EDhpf2(fc, T)
    omega = 2 * pi * fc;
    omegaD = cot(omega * T / 2);
    omegaDSq = omegaD ^ 2;
    
    b = [omegaDSq; -2 * omegaDSq; omegaDSq];
    a = [1 + sqrt(2) * omegaD + omegaDSq; 2 - 2 * omegaDSq; 1 - sqrt(2) * omegaD + omegaDSq];

%     K = tan(2 * pi * fc * T / 2);
%     K2 = K ^ 2;
%     b = [1; -2; 1];
%     a = [1 + sqrt(2) * K + K2; 2 * (K2 - 1); (1 - K * sqrt(2) + K2)];

    b = b / a(1);
    a = a / a(1);
end