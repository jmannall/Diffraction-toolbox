%% Implement delay line

function output = DelayLine(audio, pathLength, validPath, c, samplesPerUpdate, fs)

    [delay, fracDelay, amplitude] = CalculateDelay(pathLength, c, fs);

    [buffer, read, write, window, overlap, numBuffers, inputBuffer, output] = InitialiseBuffers(delay, samplesPerUpdate, audio, pathLength);
    
    disp('Process delay line')
    for k = 1:numBuffers
        [read, write] = UpdateReadWrite(read, write, overlap, delay, k);
        for i = 1:samplesPerUpdate
            [idx, read, write, inputBuffer, input] = ProcessSampleData(samplesPerUpdate, overlap, read, write, inputBuffer, amplitude, validPath, buffer, k, i);
            [inputBuffer, output, buffer] = ProcessSampleOutput(audio, buffer, input, inputBuffer, fracDelay, window, write, idx, k, i, output);
        end
    end
end