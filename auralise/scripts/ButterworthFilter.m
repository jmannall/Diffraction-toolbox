test = ButterworthFilterTesting(1000, 48000);

function output = ButterworthFilterTesting(input, fc, fs)
    K = tan(2 * pi * fc / (2 * fs));
    K2 = K ^ 2;
    a1 = 1 + sqrt(2) * K + K2;
    g = K2 / a1;
    b = [1 2 1];
    a = [1 2 * (K2 - 1) / a1 (1 - sqrt(K) + K2) / a1];
end