%% Implement delay line

function output = DelayLineBiquad(audio, pathLength, validPath, b, a, numBiquads, c, samplesPerUpdate, fs)

    [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs);
    amplitude = 2 * amplitude;

    [buffer, read, write, window, overlap, numBuffers, inputBuffer, output] = InitialiseBuffers(delay, samplesPerUpdate, audio, pathLength);
    [inputBufferIIR, outputBufferIIR] = InitialiseIIRBuffers(numBiquads);     % IIR filter input/output buffers
    
    disp('Process IIR filter and delay line')
    for k = 1:numBuffers
        [read, write] = UpdateReadWrite(read, write, overlap, delay, k);
        for i = 1:samplesPerUpdate
            [idx, read, write, inputBuffer, input] = ProcessSampleData(samplesPerUpdate, overlap, read, write, inputBuffer, amplitude, validPath, buffer, k, i);
            [inputBufferIIR, outputBufferIIR, input] = ProcessIIRFilter(input, b, a, inputBufferIIR, outputBufferIIR, numBiquads, k);
            [inputBuffer, output, buffer] = ProcessSampleOutput(audio, buffer, input, inputBuffer, fracDelay, window, write, idx, k, i, output);
        end
    end
end