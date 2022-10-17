%% Implement delay line with a biquad filter structure

function output = DelayLineBiquad(audio, pathLength, windowLength, validPath, b, a, numBiquads, c, fs)

    [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs);
    % amplitude = 2 * amplitude;    % Added as using a neural network
    % trained of radii of 1 so the overall amplitude 6dB lower than the BTM
    % comparison with radii 0.5 (total pathlength of 1)

    [buffer, read, write, window, overlap, numBuffers, inputBuffer, output] = InitialiseBuffers(delay, windowLength, audio, pathLength);
    [inputBufferIIR, outputBufferIIR] = InitialiseIIRBuffers(numBiquads);     % IIR filter input/output buffers
    
    disp('Process IIR filter and delay line')
    for k = 1:numBuffers
        [read, write] = UpdateReadWrite(read, write, overlap, delay, k);
        for i = 1:windowLength
            [idx, read, write, inputBuffer, input] = ProcessSampleData(windowLength, overlap, read, write, inputBuffer, amplitude, validPath, buffer, k, i);
            [inputBufferIIR, outputBufferIIR, input] = ProcessIIRFilter(input, b(:,:,k), a(:,:,k), inputBufferIIR, outputBufferIIR, numBiquads, k);
            [inputBuffer, output, buffer] = ProcessSampleOutput(audio, buffer, input, inputBuffer, fracDelay, window, write, idx, k, i, output);
        end
    end
end