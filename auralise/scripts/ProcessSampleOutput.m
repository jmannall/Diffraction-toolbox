function [inputBuffer, output, buffer] = ProcessSampleOutput(audio, buffer, input, inputBuffer, fracDelay, window, write, idx, k, i, output)
    
    inputBuffer(1) = input;
    %input = inputBuffer(1) * (1 - fracDelay(k)) + inputBuffer(2) * fracDelay(k); % A(1-r)*d(n-n0)+A*r*d(n-n0-1);
    output(idx) = output(idx) + window(i) * (1 - fracDelay(k)) * input;
    output(idx + 1) = output(idx + 1) + window(i + 1) * fracDelay(k) * input;
    buffer(write) = audio(idx);
end