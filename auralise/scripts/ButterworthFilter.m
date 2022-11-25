function output = ButterworthFilter(input, fc, fs, inputBufferIIR, outputBufferIIR, numBiquads, k)
    K = tan(2 * pi * fc / (2 * fs));
    K2 = K ^ 2;
    a1 = 1 + sqrt(2) * K + K2;
    g = K2 / a1;
    b = [1 2 1];
    a = [1 2 * (K2 - 1) / a1 (1 - sqrt(K) + K2) / a1];

    b = g .* [b; b];
    a = [a; a];
    numBiquads = 2;

    for k = 1:numBuffers
        [inputBufferIIR, outputBufferIIR, input] = ProcessIIRFilter(input, b, a, inputBufferIIR, outputBufferIIR, numBiquads, k);
    end
end